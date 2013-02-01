# Copyright (C) 2013 DNAnexus, Inc.
#
# This file is part of gatk_pipeline (DNAnexus platform app).
#
# (The MIT Expat License)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator
import time
import string

import dxpy
import subprocess, logging
import os, sys, re, math, operator, time
from multiprocessing import Pool, cpu_count

def main():
    ## RUN DEDUPLICATE
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar'

    #mappingsTable = dxpy.open_dxgtable(job['input']['mappings'][0]['$dnanexus_link'])
    if 'output_name' in job['input']:
        outputName = job['input']['output_name']
    else:
        outputName = ''
    recalibratedTable = createNewMappingsTable(job['input']['mappings'], outputName)
    #print "Mappings Table: " + mappingsTable.get_id()
    print "Recalibrated Table: " + recalibratedTable.get_id()

    mappingsTable = dxpy.DXGTable(job['input']['mappings'][0]['$dnanexus_link'])
    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise dxpy.AppError("The original reference genome must be attached as a detail")

    reads = 0
    for x in job['input']['mappings']:
        table = dxpy.DXGTable(x)
        reads += int(table.describe()['length'])
    chunks = int(reads/job['input']['reads_per_job'])+1

    gatkCommand = buildCommand(job)
    variantsTable = buildVariantsTable(job, mappingsTable)


    reduceInput = {}
    #Split the genome into chunks to parallelize
    commandList = splitGenomeLengthChromosome(originalContigSet, chunks)
    chunks = len(commandList)
    excludeInterchromosome = (chunks > 1)
    markDuplicatesJobs = []

    #This is a Picard Mark Duplicates job run only on interchromosomal mappings in the case that the genome is split into regions
    #This is necessary because Mark Duplicates has to look at both mates in a read pair, so interchromosomal mappings must go together
    reduceInterchromosomeInput = {}
    bamFiles = []
    if chunks > 1 and job['input']['deduplicate_interchromosomal_pairs']:
        for i in xrange(-1, chunks):
            bamFiles.append(dxpy.new_dxfile().get_id())
            mapInterchromosomeInput = {
            'mappings_tables': job['input']['mappings'],
            'interval': commandList[i],
            'job_number' : i,
            'separate_read_groups' : job['input']['separate_read_groups']
            }
            interchromosomeJobId = dxpy.new_dxjob(fn_input=mapInterchromosomeInput, fn_name="mapInterchromosome").get_id()
            reduceInterchromosomeInput["mapJob"+str(i)] = {'job': interchromosomeJobId, 'field': 'file_id'}
        #interchromosomeJobField = { 'job': interchromosomeJobId, 'field': 'bam'}

        #Make a reduce job for the interchromosome component
        reduceInterchromosomeInput["file_list"] =  bamFiles
        reduceInterchromosomeInput["interval"] = commandList
        reduceInterchromosomeInput["discard_duplicates"] = job['input']['discard_duplicates']
        reduceJobId = dxpy.new_dxjob(fn_input=reduceInterchromosomeInput, fn_name="reduceInterchromosome").get_id()
        deduplicateInterchromosome = True
    else:
        interchromosomeJobField = ''
        deduplicateInterchromosome = False

    #This runs the Picard Mark Duplicates program to deduplicate the reads
    reduceInput = {}
    variantCallingCoordinatorInput = {}
    for i in range(len(commandList)):
        mapBestPracticesInput = {
            'mappings_tables': job['input']['mappings'],
            'recalibrated_table_id': recalibratedTable.get_id(),
            'file_list': bamFiles,
            'interval': commandList[i],
            'job_number' : i,
            'reference': job['input']['reference']['$dnanexus_link'],
            'dbsnp': job['input']['dbsnp'],
            'separate_read_groups' : job['input']['separate_read_groups'],
            'discard_duplicates': job['input']['discard_duplicates'],
            'parent_input': job['input'],
            'deduplicate_interchromosome': deduplicateInterchromosome,
            'gatk_command': gatkCommand,
            'compress_reference': job['input']['compress_reference'],
            'infer_no_call': job['input']['infer_no_call'],
            'compress_no_call': job['input']['compress_no_call'],
            'variants_table': variantsTable.get_id()
        }
        if 'known_indels' in job['input']:
            mapBestPracticesInput['known_indels'] = job['input']['known_indels']


        mapJobId = dxpy.new_dxjob(fn_input=mapBestPracticesInput, fn_name="mapBestPractices").get_id()
        variantCallingCoordinatorInput['best_practices_job'+str(i)] = {'job': mapJobId, 'field': 'gatk_job'}
        reduceInput["mapJob" + str(i)] = {'job': mapJobId, 'field': 'ok'}
    variantCallingCoordinatorJob = dxpy.new_dxjob(fn_input=variantCallingCoordinatorInput, fn_name="variantCallingCoordinator").get_id()
    reduceInput["recalibrated_table"] = recalibratedTable.get_id()
    reduceInput["variants_table"] = variantsTable.get_id()
    reduceInput["gatk_completion"] = {'job': variantCallingCoordinatorJob, 'field': 'gatk_jobs'}


    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceBestPractices").get_id()
    job['output'] = {'recalibrated_mappings': {'job': reduceJobId, 'field': 'recalibrated_table'}, 'variants': {'job': reduceJobId, 'field': 'variants_table'}}

def variantCallingCoordinator():
    job['output']['gatk_jobs'] = []
    for i in range(len(job['input'])):
        job['output']['gatk_jobs'].append({'job': job['input']['best_practices_job'+str(i)], 'field': 'ok'})


def reduceBestPractices():
    startTime = time.time()
    recalibratedTable = dxpy.DXGTable(job['input']['recalibrated_table'])
    recalibratedTable.close()
    print "Table closing completed in " + str(int((time.time()-startTime)/60)) + " minutes"
    variantsTable = dxpy.DXGTable(job['input']['variants_table'])
    variantsTable.close()
    job['output']['recalibrated_table'] = dxpy.dxlink(recalibratedTable.get_id())
    job['output']['variants_table'] = dxpy.dxlink(variantsTable.get_id())

def writeUnmappedReads(mappingsTable, dedupTable):
    colNames = dedupTable.get_col_names()
    col = {}
    for i in range(len(colNames)):
        col[colNames[i]] = i
    for row in mappingsTable.iterate_rows():
        entry = []
        if row[col["chr"]] != '':
            break
        if row[col["chr"]] == row[col["chr2"]]:
            dedupTable.add_rows([entry])


def buildVariantsTable(job, mappingsTable):
    if job['input']['output_mode'] == "EMIT_VARIANTS_ONLY":
        job['input']['infer_no_call'] = False

    variants_schema = [
        {"name": "chr", "type": "string"},
        {"name": "lo", "type": "int32"},
        {"name": "hi", "type": "int32"},
        {"name": "ref", "type": "string"},
        {"name": "alt", "type": "string"},
        {"name": "qual", "type": "double"},
        {"name": "ids", "type": "string"}
         ]

    elevatedTags = ['format_GT', 'format_DP', 'format_AD']
    headerInfo = extractHeader("/tmp/header.txt", elevatedTags)
    description = {}
    samples = []

    print headerInfo
    print headerInfo['tags']['format']

    indices = [dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')]

    formats = {}
    infos = {}
    filters = {}

    for k, v in headerInfo['tags']['info'].iteritems():
        variants_schema.append({"name": "info_"+k, "type":translateTagTypeToColumnType(v)})
        description[k] = {'name' : k, 'description' : v['description'], 'type' : v['type'], 'number' : v['number']}

    numSamples = 1
    #For each sample, write the sample-specific columns
    for i in range(numSamples):
      variants_schema.extend([
        {"name": "genotype_"+str(i), "type": "string"},
        {"name": "phasing_"+str(i), "type": "string"},
        {"name": "type_"+str(i), "type": "string"},
        {"name": "variation_qual_"+str(i), "type": "double"},
        {"name": "genotype_qual_"+str(i), "type": "double"},
        {"name": "coverage_"+str(i), "type": "string"},
        {"name": "total_coverage_"+str(i), "type": "int32"}
      ])
      samples.append("Sample_0")
      for k, v in headerInfo['tags']['format'].iteritems():
        if "format_"+k not in elevatedTags:
          variants_schema.append({"name": "format_"+k+"_"+str(i), "type":translateTagTypeToColumnType(v)})

    variantsTable = dxpy.new_dxgtable(variants_schema, indices=[dxpy.DXGTable.genomic_range_index("chr", "lo", "hi", "gri")])
    tableId = variantsTable.get_id()
    variantsTable = dxpy.open_dxgtable(tableId)
    variantsTable.add_types(["Variants", "gri"])

    details = {'samples':samples, 'original_contigset':job['input']['reference'], 'original_mappings':job['input']['mappings'], 'formats':headerInfo['tags']['format'], 'infos':headerInfo['tags']['info']}
    #if headerInfo.get('filters') != {}:
    #  details['filters'] = headerInfo['filters']
    variantsTable.set_details(details)

    if 'output_name' in job['input']:
        variantsTable.rename(job['input']['output_name'])
    elif (job['input']['genotype_likelihood_model'] == "SNP"):
        variantsTable.rename(mappingsTable.describe()['name'] + " SNP calls by GATK")
    elif (job['input']['genotype_likelihood_model'] == "INDEL"):
        variantsTable.rename(mappingsTable.describe()['name'] + " indel calls by GATK")
    elif (job['input']['genotype_likelihood_model'] == "BOTH"):
        variantsTable.rename(mappingsTable.describe()['name'] + " SNP and indel calls by GATK")
    else:
        variantsTable.rename(mappingsTable.describe()['name'] + " variant calls by GATK")

    return variantsTable

def mapInterchromosome():
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/AddOrReplaceReadGroups.jar'

    regionFile = open("regions.txt", 'w')
    regionFile.write(job['input']['interval'])
    regionFile.close()

    jobNumber = job['input']['job_number']
    if jobNumber == -1:
        for i in range(len(job['input']['mappings_tables'])):
            mappingsTable = dxpy.DXGTable(job['input']['mappings_tables'][i]['$dnanexus_link']).get_id()
            command = "pypy /usr/bin/dx_mappings_to_sam2 %s --output input.%d.sam --id_as_name --only_unmapped --read_group_platform illumina --id_prepend %d_" % (mappingsTable, i, i)
            print command
            subprocess.check_call(command, shell=True)
    else:
        for i in range(len(job['input']['mappings_tables'])):
            mappingsTable = dxpy.DXGTable(job['input']['mappings_tables'][i]['$dnanexus_link']).get_id()
            command = "pypy /usr/bin/dx_mappings_to_sam2 %s --output input.%d.sam --id_as_name --region_index_offset -1 --region_file regions.txt --read_group_platform illumina --only_interchromosomal_mate --id_prepend %d_" % (mappingsTable, i, i)
            print command
            subprocess.check_call(command, shell=True)


    readsPresent = False
    if len(job['input']['mappings_tables']) == 1:
        if checkSamContainsRead("input.0.sam"):
            readsPresent = True
            subprocess.check_call("mv input.0.sam output.sam", shell=True)
    else:
        command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=output.sam USE_THREADING=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
        for i in range(len(job['input']['mappings_tables'])):
            if checkSamContainsRead("input."+str(i)+".sam"):
                command += " INPUT=input."+str(i)+".sam"
                readsPresent = True
        subprocess.check_call(command, shell=True)
    if readsPresent:
        subprocess.check_call("samtools view -bS output.sam > output.bam", shell=True)
        fileId = dxpy.upload_local_file("output.bam").get_id()
        job['output']['file_id'] = fileId
    else:
        job['output']['file_id'] = ''

def reduceInterchromosome():
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/AddOrReplaceReadGroups.jar'
    i = -1
    command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=merge.bam USE_THREADING=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
    noFiles = True
    while "mapJob"+str(i) in job["input"]:
        if job["input"]["mapJob"+str(i)] != '':
            command += " INPUT=input."+str(i+1)+".bam"
            dxpy.DXFile(job["input"]["mapJob"+str(i)]).wait_on_close()
            dxpy.download_dxfile(job["input"]["mapJob"+str(i)], "input."+str(i+1)+".bam")
            noFiles = False
        i+=1

    if noFiles:
        for x in job['input']['file_list']:
            dxpy.DXFile(x).close()
    else:
        subprocess.check_call(command, shell=True)
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=merge.bam O=dedup.bam METRICS_FILE=metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=%s" % job['input']['discard_duplicates'], shell=True)
        #subprocess.check_call("samtools view -bS dedup.sam > dedup.bam", shell=True)
        subprocess.check_call("samtools index dedup.bam", shell=True)
        for i in range(len(job['input']['interval'])):
            try:
                startTime = time.time()
                chromosomeFile = open("chromosome.bam", 'w')
                chromosomeFile.close()
                fileId = dxpy.DXFile(job['input']['file_list'][i])
                command = "samtools view -b dedup.bam " + job['input']['interval'][i].replace(" -L ", " ") + " > chromosome.bam"
                print command
                subprocess.check_call(command, shell=True)
                if os.stat('chromosome.bam').st_size > 0:
                    dxpy.upload_local_file("chromosome.bam", use_existing_dxfile=fileId)
                    print "Result upload completed in " + str(int((time.time()-startTime)/60)) + " minutes"
                else:
                    print "This interval was empty of reads"
                    dxpy.DXFile(job['input']['file_list'][i]).close()
            except:
                print "Could not write this interval"
                dxpy.DXFile(job['input']['file_list'][i]).close()
    job['output']['ok'] = True
    return

def mapBestPractices():
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/AddOrReplaceReadGroups.jar'

    jobNumber = job['input']['job_number']

    regionFile = open("regions.txt", 'w')
    print job['input']['interval']
    regionFile.write(job['input']['interval'])
    regionFile.close()

    readGroups = 0
    print "Converting Table to SAM"
    for i in range(len(job['input']['mappings_tables'])):
        mappingsTable = dxpy.DXGTable(job['input']['mappings_tables'][i]['$dnanexus_link']).get_id()
        command = "pypy /usr/bin/dx_mappings_to_sam2 %s --output input.%d.sam --region_index_offset -1 --id_as_name --region_file regions.txt --write_row_id --read_group_platform illumina --id_prepend %d_" % (mappingsTable, i, i)
        if job['input']['deduplicate_interchromosome']:
            command += " --no_interchromosomal_mate"
        if job['input']['separate_read_groups']:
            command += " --add_to_read_group " + str(readGroups)
            readGroups += len(dxpy.DXGTable(job['input']['mappings_tables'][i]['$dnanexus_link']).get_details()['read_groups'])
        print command
        startTime = time.time()
        subprocess.check_call(command, shell=True)
        print "Download mappings completed in " + str(int((time.time()-startTime)/60)) + " minutes"

    readsPresent = False

    if len(job['input']['mappings_tables']) == 1:
        if checkSamContainsRead("input.0.sam"):
            readsPresent = True
            subprocess.check_call("mv input.0.sam input.sam", shell=True)
    else:
        command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=input.sam USE_THREADING=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
        for i in range(len(job['input']['mappings_tables'])):
            if checkSamContainsRead("input."+str(i)+".sam"):
                command += " INPUT=input."+str(i)+".sam"
                readsPresent = True

    if readsPresent:
        subprocess.check_call(command, shell=True)
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=input.sam O=dedup.sam METRICS_FILE=metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=%s" % job["input"]["discard_duplicates"], shell=True)
        startTime = time.time()
        subprocess.check_call("samtools view -bS dedup.sam > dedup.bam", shell=True)
        print "Conversion to BAM completed in " + str(int((time.time()-startTime)/60)) + " minutes"
    else:
        job['output']['gatk_job'] = True
        job['output']['ok'] = True
        return

    ##THIS SEGMENT PROCEEDS AFTER MARKDUPLICATES COMPLETES
    sleepTime = 10
    sleepCounter = 0
    startTime = time.time()
    if jobNumber >  -1 and job['input']['deduplicate_interchromosome']:
        if job['input']['file_list'][jobNumber] != '':
            interchromosomeBamId = job['input']['file_list'][jobNumber]
            interchromosomeBam = dxpy.DXFile(interchromosomeBamId)
            while 1:
                if interchromosomeBam.describe()['state'] != 'closed':
                    time.sleep(sleepTime)
                    sleepCounter += sleepTime
                else:
                    print "Waited " + str(sleepTime/60) + " minutes for interchromosome job"
                    dxpy.download_dxfile(interchromosomeBamId, "interchromosomeBam.bam")
                    print "Interchromosome BAM: " + interchromosomeBam.get_id()
                    print "Download interchromosome BAM in " + str(int((time.time()-startTime)/60)) + " minutes"
                    break

    print "dedup file:" + dxpy.upload_local_file("dedup.bam").get_id()

    if job['input']['file_list'] == []:
        subprocess.check_call("mv dedup.bam input.bam", shell=True)
    elif job['input']['file_list'][jobNumber] != '':
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.MergeSamFiles SORT_ORDER=coordinate USE_THREADING=true INPUT=dedup.bam INPUT=interchromosomeBam.bam OUTPUT=input.bam VALIDATION_STRINGENCY=SILENT", shell=True)
    else:
        subprocess.check_call("mv dedup.bam input.bam", shell=True)
    startTime = time.time()
    subprocess.check_call("samtools index input.bam", shell=True)
    print "Index BAM completed in " + str(int((time.time()-startTime)/60)) + " minutes"

    #Download the Reference Genome
    print "Converting Contigset to Fasta"
    subprocess.check_call("contigset2fasta %s ref.fa" % (job['input']['reference']), shell=True)

    #RealignerTargetCreator
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T RealignerTargetCreator -R ref.fa -I input.bam -o indels.intervals "
    command += job['input']['interval']
    knownIndels = ''

    #Download the known indels
    knownCommand = ''
    if 'known_indels' in job['input']:
        for i in range(len(job['input']['known_indels'])):
            dxpy.download_dxfile(job['input']['known_indels'][i], "indels"+str(i)+".vcf.gz")
            knownFileName = "indels"+str(i)+".vcf.gz"
            try:
                p = subprocess.Popen("tabix -f -p vcf " + knownFileName, stderr=subprocess.PIPE, shell=True)
                if '[tabix] was bgzip' in p.communicate()[1]:
                    subprocess.check_call("gzip -d " + knownFileName, shell=True)
                    knownFileName = "indels"+str(i)+".vcf.gz"
            except subprocess.CalledProcessError:
                raise dxpy.AppError("An error occurred decompressing dbsnp. Expected dbsnp as a gziped file.")
            knownCommand += " -known " + knownFileName
        command += knownIndels

    #Find chromosomes
    regionChromosomes = []
    for x in re.findall("-L ([^:]*):\d+-\d+", job['input']['interval']):
        regionChromosomes.append(x)

    #Add options for RealignerTargetCreator
    if job['input']['parent_input']['window_size'] != 10:
        command += "--windowSize " + str(job['input']['parent_input']['window_size'])
    if job['input']['parent_input']['max_interval_size'] != 500:
        command += " --maxIntervalSize "  + str(job['input']['parent_input']['max_interval_size'])
    if job['input']['parent_input']['min_reads_locus'] != 4:
        command += " --minReadsAtLocus " + str(job['input']['parent_input']['min_reads_locus'])
    if job['input']['parent_input']['mismatch_fraction'] != 0.0:
        command += " --mismatchFraction " + str(job['input']['parent_input']['mismatch_fraction'])

    print command
    subprocess.check_call(command, shell=True)

    #Run the IndelRealigner
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T IndelRealigner -R ref.fa -I input.bam -targetIntervals indels.intervals -o realigned.bam"
    command += job['input']['interval']
    command += knownCommand
    if "consensus_model" in job['input']['parent_input']:
        if job['input']['parent_input']['consensus_model'] != "":
            if job['input']['parent_input']['consensus_model'] == "USE_READS" or job['input']['parent_input']['consensus_model'] == "KNOWNS_ONLY" or job['input']['parent_input']['consensus_model'] == "USE_SW":
                command += " --consensusDeterminationModel " + job['input']['parent_input']['consensus_model']
            else:
                raise dxpy.AppError("The option \"Consensus Determination Model\" must be either blank or one of [\"USE_READS\", \"KNWONS_ONLY\", or \"USE_SW\"], found " + job['input']['parent_input']['consensus_model'] + " instead.")
    if job['input']['parent_input']['lod_threshold'] != 5.0:
        command += " --LODThresholdForCleaning " + str(job['input']['parent_input']['lod_threshold'])
    if job['input']['parent_input']['entropy_threshold'] != 0.15:
        command += " --entropyThreshold " + str(job['input']['parent_input']['entropy_threshold'])
    if job['input']['parent_input']['max_consensuses'] != 30:
        command += " --maxConsensuses " + str(job['input']['parent_input']['max_consensuses'])
    if job['input']['parent_input']['max_insert_size_movement'] != 3000:
        command += " --maxIsizeForMovement " + str(job['input']['parent_input']['max_insert_size_movement'])
    if job['input']['parent_input']['max_position_move'] != 200:
        command += " --maxPositionalMoveAllowed " + str(job['input']['parent_input']['max_position_move'])
    if job['input']['parent_input']['max_reads_consensus'] != 120:
        command += " --maxReadsForConsensus " + str(job['input']['parent_input']['maxReadsForRealignment'])
    if job['input']['parent_input']['max_reads_realignment'] != 20000:
        command += " --maxReadsForRealignment " + str(job['input']['parent_input']['max_reads_realignment'])

    print command
    subprocess.check_call(command, shell=True)

    #Download dbsnp
    startTime = time.time()
    dxpy.download_dxfile(job['input']['dbsnp'], "dbsnp.vcf.gz")
    print "Download dbsnp completed in " + str(int((time.time()-startTime)/60)) + " minutes"

    dbsnpFileName = 'dbsnp.vcf.gz'
    try:
        p = subprocess.Popen("tabix -f -p vcf dbsnp.vcf.gz", stderr=subprocess.PIPE, shell=True)
        if '[tabix] was bgzip' in p.communicate()[1]:
            subprocess.check_call("gzip -d dbsnp.vcf.gz", shell=True)
            dbsnpFileName = 'dbsnp.vcf'
    except subprocess.CalledProcessError:
        raise dxpy.AppError("An error occurred decompressing dbsnp. Expected dbsnp as a gziped file.")

    #Count Covariates
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T CountCovariates -R ref.fa -recalFile recalibration.csv -I realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate --standard_covs"
    command += " -knownSites " + dbsnpFileName
    command += job['input']['interval']
    command += " --num_threads " + str(cpu_count())

    if "solid_recalibration_mode" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_recalibration_mode'] != "":
            if job['input']['parent_input']['solid_recalibration_mode'] == "THROW_EXCEPTION" or job['input']['parent_input']['solid_recalibration_mode'] == "LEAVE_READ_UNRECALIBRATED" or job['input']['parent_input']['solid_recalibration_mode'] == "PURGE_READ":
                command += " --solid_recal_mode " + job['input']['parent_input']['solid_recalibration_mode']
            else:
                raise dxpy.AppError("The option \"SOLiD Recalibration Mode\" must be either blank or one of [\"THROW_EXCEPTION\", \"LEAVE_READ_UNRECALIBRATED\", or \"PURGE_READ\"], found " + job['input']['parent_input']['solid_recalibration_mode'] + " instead.")
    if "solid_nocall_strategy" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_nocall_strategy'] != "":
            if job['input']['parent_input']['solid_nocall_strategy'] == "THROW_EXCEPTION" or job['input']['parent_input']['solid_nocall_strategy'] == "LEAVE_READ_UNRECALIBRATED" or job['input']['parent_input']['solid_nocall_strategy'] == "PURGE_READ":
                command += " --solid_nocall_strategy " + job['input']['parent_input']['solid_nocall_strategy']
            else:
                raise dxpy.AppError("The option \"SOLiD No-call Strategy\" must be either blank or one of [\"THROW_EXCEPTION\", \"LEAVE_READ_UNRECALIBRATED\", or \"PURGE_READ\"], found " + job['input']['parent_input']['solid_nocall_strategy'] + " instead.")
    if 'context_size' in job['input']['parent_input']:
        command += " --context_size " + str(job['input']['parent_input']['context_size'])
    if 'nback' in job['input']['parent_input']:
        command += " --homopolymer_nback " + str(job['input']['parent_input']['nback'])
    if job['input']['parent_input']['cycle_covariate']:
        command += " -cov CycleCovariate"
    if job['input']['parent_input']['dinuc_covariate']:
        command += " -cov DinucCovariate"
    if job['input']['parent_input']['primer_round_covariate']:
        command += " -cov PrimerRoundCovariate"
    if job['input']['parent_input']['mapping_quality_covariate']:
        command += " -cov MappingQualityCovariate"
    if job['input']['parent_input']['homopolymer_covariate']:
        command += " -cov HomopolymerCovariate"
    if job['input']['parent_input']['gc_content_covariate']:
        command += " -cov GCContentCovariate"
    if job['input']['parent_input']['position_covariate']:
        command += " -cov PositionCovariate"
    if job['input']['parent_input']['minimum_nqs_covariate']:
        command += " -cov MinimumNQSCovariate"
    if job['input']['parent_input']['context_covariate']:
        command += " -cov ContextCovariate"

    print command
    subprocess.check_call(command, shell=True)

    #Table Recalibration
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T TableRecalibration -R ref.fa -recalFile recalibration.csv -I realigned.bam -o recalibrated.bam --doNotWriteOriginalQuals"
    command += job['input']['interval']
    if "solid_recalibration_mode" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_recalibration_mode'] != "":
            command += " --solid_recal_mode " + job['input']['parent_input']['solid_recalibration_mode']
    if "solid_nocall_strategy" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_nocall_strategy'] != "":
            command += " --solid_nocall_strategy " + job['input']['parent_input']['solid_nocall_strategy']
    if 'context_size' in job['input']['parent_input']:
        command += " --context_size " + str(job['input']['parent_input']['context_size'])
    if 'nback' in job['input']['parent_input']:
        command += " --homopolymer_nback " + str(job['input']['parent_input']['nback'])
    if job['input']['parent_input']['preserve_qscore'] != 5:
        command += " --preserve_qscores_less_than " + str(job['input']['parent_input']['preserve_qscore'])
    if 'smoothing' in job['input']['parent_input']:
        command += " --smoothing " + str(job['input']['parent_input']['smoothing'])
    if 'max_quality' in job['input']['parent_input']:
        command += " --max_quality_score " + str(job['input']['parent_input']['max_quality'])
    print command
    subprocess.check_call(command, shell=True)

    subprocess.check_call("samtools view -h -o recalibrated.sam recalibrated.bam", shell=True)

    result = dxpy.upload_local_file("recalibrated.bam")
    job['output']['recalibrated_bam'] = dxpy.dxlink(result.get_id())
    print "Recalibrated file: " + result.get_id()
    referenceFile = dxpy.upload_local_file("ref.fa")

    #Spin off GATK variant caller job
    mapGatkInput = {
            'mappings_file': result.get_id(),
            'reference_file': referenceFile.get_id(),
            'interval': job['input']['interval'],
            'tableId': job['input']['variants_table'],
            'command': job['input']['gatk_command'],
            'compress_reference': job['input']['compress_reference'],
            'infer_no_call': job['input']['infer_no_call'],
            'compress_no_call': job['input']['compress_no_call'],
            'part_number': job['input']['job_number']
        }
    # Run a "map" job for each chunk
    gatkJob = dxpy.new_dxjob(fn_input=mapGatkInput, fn_name="mapGatk").get_id()

    recalibratedTable = dxpy.DXGTable(job['input']['recalibrated_table_id'])

    default = {}

    recalibratedColNames = recalibratedTable.get_col_names()
    recalibratedCol = {}
    for i in range(len(recalibratedColNames)):
        recalibratedCol[recalibratedColNames[i]] = i

    for x in recalibratedTable.describe()['columns']:
        if "int" in x["type"]:
            default[x["name"]] = sys.maxint
        elif x["type"] == "float":
            default[x["name"]] = float(sys.maxint)
        elif x["type"] == "boolean":
            default[x["name"]] = False
        else:
            default[x["name"]] = ""

    print "Writing mate pair information for lookup"
    startTime = time.time()
    if recalibratedCol.get("chr2") != None:
        mateLocations = {}
        for line in open("recalibrated.sam", 'r'):
            if line[0] != "@":
                tabSplit = line.split("\t")
                chr = tabSplit[2]
                lo = int(tabSplit[3])-1
                templateId = tabSplit[0]
                cigar = re.split('(\d+)', tabSplit[5])
                alignLength = 0
                for p in range(len(cigar)):
                    c = cigar[p]
                    if c == 'M' or c == 'D' or c == 'N' or c == 'X' or c == 'P' or c == '=':
                        alignLength += int(cigar[p-1])
                hi = lo + alignLength

                recalibrationTags = re.findall("zd:Z:([^\t]*)[\t\n]", line)[0].split("##&##")
                reportedLo = int(recalibrationTags[1])
                reportedHi = int(recalibrationTags[2])
                if lo != reportedLo or hi != reportedHi:
                    if int(tabSplit[1]) & 0x1 & 0x40:
                        mateLocations[templateId] = {0: {"lo":lo, "hi":hi, "chr":chr}}
                    elif int(tabSplit[1]) & 0x1 & 0x80:
                        mateLocations[templateId] = {1: {"lo":lo, "hi":hi, "chr":chr}}
        print str(len(mateLocations)) + " Interchromosomal reads changed lo or hi"
    print "Interchromosome changes to lo and hi recorded in " + str(int((time.time()-startTime)/60)) + " minutes"

    complement_table = string.maketrans("ATGCatgc", "TACGtacg")
    rowsWritten = 0

    startTime = time.time()
    for line in open("recalibrated.sam", 'r'):
        if line[0] != "@":
            tabSplit = line.split("\t")
            templateId = int(tabSplit[0].split("_")[1])
            flag = int(tabSplit[1])
            chr = tabSplit[2]
            lo = int(tabSplit[3])-1
            qual = tabSplit[10]
            alignLength = 0
            mapq = int(tabSplit[4])
            cigar = re.split('(\d+)', tabSplit[5])
            duplicate = (flag & 0x400 == True)
            sequence = tabSplit[9]
            recalibrationTags = re.findall("zd:Z:([^\t]*)[\t\n]", line)[0].split("##&##")

            name = recalibrationTags[0]
            readGroup = int(re.findall("RG:Z:(\d+)", line)[0])

            if flag & 0x4:
                status = "UNMAPPED"
            elif flag & 0x100:
                status = "SECONDARY"
            else:
                status = "PRIMARY"
            if flag & 0x200:
                qcFail = True
            else:
                qcFail = False
            if flag & 0x10 == 0 or status == "UNMAPPED":
                negativeStrand = False
            else:
                negativeStrand = True
                sequence = sequence.translate(complement_table)[::-1]
                qual = qual[::-1]

            for p in range(len(cigar)):
                c = cigar[p]
                if c == 'M' or c == 'D' or c == 'N' or c == 'X' or c == 'P' or c == '=':
                    alignLength += int(cigar[p-1])
            cigar = tabSplit[5]
            hi = lo + alignLength

            if recalibratedCol.get("chr2") != None:
                properPair=False
                if (flag & 0x1) and (flag & 0x2):
                    properPair = True

                reportedChr2 = recalibrationTags[3]
                reportedLo2 = int(recalibrationTags[4])
                reportedHi2 = int(recalibrationTags[5])
                status2 = recalibrationTags[6]

                if not flag & 0x1:
                    negativeStrand2 = False
                elif flag & 0x20 and status2 != "UNMAPPED":
                    negativeStrand2 = True
                else:
                    negativeStrand2 = False

                if flag & 0x1:
                    if flag & 0x40:
                        mateId = 0
                    elif flag & 0x80:
                        mateId = 1
                else:
                    mateId = -1

                try:
                    if mateId == 1:
                        chr2 = mateLocations[tabSplit[0]][0]["chr2"]
                        lo2 = mateLocations[tabSplit[0]][0]["lo"]
                        hi2 = mateLocations[tabSplit[0]][0]["hi"]
                    elif mateId == 0:
                        chr2 = mateLocations[tabSplit[0]][1]["chr2"]
                        lo2 = mateLocations[tabSplit[0]][1]["lo"]
                        hi2 = mateLocations[tabSplit[0]][1]["hi"]
                    else:
                        chr2 = reportedChr2
                        lo2 = reportedLo2
                        hi2 = reportedHi2
                except:
                    chr2 = reportedChr2
                    lo2 = reportedLo2
                    hi2 = reportedHi2
                recalibratedTable.add_rows([[sequence, name, qual, status, chr, lo, hi, negativeStrand, mapq, qcFail, duplicate, cigar, templateId, readGroup, mateId, status2, chr2, lo2, hi2, negativeStrand2, properPair]])
            else:
                recalibratedTable.add_rows([[sequence, name, qual, status, chr, lo, hi, negativeStrand, mapq, qcFail, duplicate, cigar, templateId, readGroup]])
            rowsWritten += 1
            if rowsWritten%100000 == 0:
                print "Imported " + str(rowsWritten) + " rows. Time taken: " + str(int((time.time()-startTime)/60)) + " minutes"
                recalibratedTable.flush()
    recalibratedTable.flush()

    job['output']['gatk_job'] = gatkJob
    job['output']['ok'] = True

def splitGenomeLengthChromosome(contig_set, chunks):
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

    commandList = []
    chunkSizes = []

    totalSize = sum(sizes)
    readsPerChunk = int(float(totalSize)/chunks)

    for i in range(chunks):
        commandList.append('')
        chunkSizes.append(0)

    chromosome = 0
    position = 0
    totalSizes = []
    print chunkSizes
    for x in range(chunks):
        totalSizes.append(0)
    while chromosome < len(names):
        try:
            minimum = min([x for x in chunkSizes if x > 0])
            position = chunkSizes.index(minimum)
        except:
            position = 0
        if chunkSizes[position] + sizes[chromosome] > readsPerChunk:
            try:
                position = chunkSizes.index(0)
            except:
                position = chunkSizes.index(min([x for x in chunkSizes if x > 0]))
        commandList[position] += " -L %s:%d-%d" % (names[chromosome], 1, sizes[chromosome])
        chunkSizes[position] += sizes[chromosome]
        totalSizes[position] += sizes[chromosome]
        chromosome += 1
    commandList = filter(None, commandList)
    return commandList

def checkIntervalRange(includeList, chromosome, lo, hi):
    included = False
    command = ''
    if len(includeList) == 0:
        return " -L %s:%d-%d" % (chromosome, lo, hi)
    if includeList.get(chromosome) != None:
        for x in includeList[chromosome]:
            min = lo
            max = hi
            if (lo >= x[0] and lo <= x[1]) or (hi <= x[1] and hi >= x[0]):
                if lo >= x[0] and lo <= x[1]:
                    min = lo
                elif lo <= x[0]:
                    min = x[0]
                if hi <= x[1] and hi >= x[0]:
                    max = hi
                elif hi >= x[1]:
                    max = x[1]
                command += " -L %s:%d-%d" % (chromosome, min, max)
    return command

def checkSamContainsRead(samFileName):
    for line in open(samFileName, 'r'):
        if line[0] != "@":
            return True
    return False


def createNewMappingsTable(mappingsArray, recalibratedName):

    columns = []
    tags = []
    indices = []
    types = []
    read_groups = []
    for i in range(len(mappingsArray)):
        oldTable = dxpy.DXGTable(mappingsArray[i]['$dnanexus_link'])
        if oldTable.get_details()['read_groups'] != None:
            read_groups.extend(oldTable.get_details()['read_groups'])
        for x in oldTable.describe()['columns']:
            if x not in columns:
                columns.append(x)
        for x in oldTable.describe()['indices']:
            if x not in indices:
                indices.append(x)
        for x in oldTable.describe()['tags']:
            if x not in tags:
                tags.append(x)
        for x in oldTable.describe()['types']:
            if x not in types:
                types.append(x)

    print "combined columns"
    print columns

    schema = [{"name": "sequence", "type": "string"}]
    schema.append({"name": "name", "type": "string"})
    schema.append({"name": "quality", "type": "string"})
    schema.extend([{"name": "status", "type": "string"},
                          {"name": "chr", "type": "string"},
                          {"name": "lo", "type": "int32"},
                          {"name": "hi", "type": "int32"},
                          {"name": "negative_strand", "type": "boolean"},
                          {"name": "error_probability", "type": "uint8"},
                          {"name": "qc_fail", "type": "boolean"},
                          {"name": "duplicate", "type": "boolean"},
                          {"name": "cigar", "type": "string"},
                          {"name": "template_id", "type": "int64"},
                          {"name": "read_group", "type": "int32"}])
    if {"type":"string", "name":"chr2"} in columns:
        schema.extend([{"name": "mate_id", "type": "int32"}, # TODO: int8
                              {"name": "status2", "type": "string"},
                              {"name": "chr2", "type": "string"},
                              {"name": "lo2", "type": "int32"},
                              {"name": "hi2", "type": "int32"},
                              {"name": "negative_strand2", "type": "boolean"},
                              {"name": "proper_pair", "type": "boolean"}])

    oldTable = dxpy.DXGTable(mappingsArray[0]['$dnanexus_link'])
    if recalibratedName == '':
        recalibratedName = oldTable.describe()['name'] + " Realigned and Recalibrated"
    details = oldTable.get_details()
    details['read_groups'] = read_groups
    newTable = dxpy.new_dxgtable(columns=schema, indices=indices)
    newTable.add_tags(tags)
    newTable.set_details(details)
    newTable.add_types(types)
    newTable.rename(recalibratedName)

    return newTable

def buildCommand(job):

    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -I input.bam -o output.vcf "
    if job['input']['output_mode'] != "EMIT_VARIANTS_ONLY":
        command += " -out_mode " + (job['input']['output_mode'])
    if job['input']['call_confidence'] != 30.0:
        command += " -stand_call_conf " +str(job['input']['call_confidence'])
    if job['input']['emit_confidence'] != 30.0:
        command += " -stand_emit_conf " +str(job['input']['emit_confidence'])
    if job['input']['pcr_error_rate'] != 0.0001:
        command += " -pcr_error " +str(job['input']['pcr_error_rate'])
    if job['input']['heterozygosity'] != 0.001:
        command += " -hets " + str(job['input']['heterozygosity'])
    if job['input']['indel_heterozygosity'] != 0.000125:
        command += " -indelHeterozygosity " + str(job['input']['indel_heterozygosity'])
    if job['input']['genotype_likelihood_model'] != "SNP":
        command += " -glm " + job['input']['genotype_likelihood_model']
    if job['input']['minimum_base_quality'] != 17:
        command += " -mbq " + str(job['input']['minimum_base_quality'])
    if job['input']['max_alternate_alleles'] != 3:
        command += " -maxAlleles " + str(job['input']['max_alternate_alleles'])
    if job['input']['max_deletion_fraction'] != 0.05:
        command += " -deletions " + str(job['input']['max_deletion_fraction'])
    if job['input']['min_indel_count'] != 5:
        command += " -minIndelCnt " + str(job['input']['min_indel_count'])
    if job['input']['non_reference_probability_model'] != "EXACT":
        if job['input']['non_reference_probability_model'] != "GRID_SEARCH":
            raise dxpy.AppError("Option \"Probability Model\" must be either \"EXACT\" or \"GRID_SEARCH\". Found " + job['input']['non_reference_probability_model'] + " instead")
        command += " -pnrm " + str(job['input']['non_reference_probability_model'])

    command += " --num_threads " + str(cpu_count())
    command += " -L regions.interval_list"

    if job['input']['downsample_to_coverage'] != 250:
        command += " -dcov " + str(job['input']['downsample_to_coverage'])
    elif job['input']['downsample_to_fraction'] != 1.0:
        command += " -dfrac " + str(job['input']['downsample_to_fraction'])

    if job['input']['nondeterministic']:
        command += " -ndrs "

    if job['input']['calculate_BAQ'] != "OFF":
        if job['input']['calculate_BAQ'] != "CALCULATE_AS_NECESSARY" and job['input']['calculate_BAQ'] != "RECALCULATE":
            raise dxpy.AppError("Option \"Calculate BAQ\" must be either \"OFF\" or or \"CALCULATE_AS_NECESSARY\" \"RECALCULATE\". Found " + job['input']['calculate_BAQ'] + " instead")
        command += "-baq " + job['input']['calculate_BAQ']
        if job['input']['BAQ_gap_open_penalty'] != 40.0:
            command += "-baqGOP " + str(job['input']['BAQ_gap_open_penalty'])
    if job['input']['no_output_SLOD']:
        command += "-nosl"

    return command

def extractHeader(vcfFileName, elevatedTags):
    result = {'columns': '', 'tags' : {'format' : {}, 'info' : {} }, 'filters' : {}}
    for line in open(vcfFileName):
        tag = re.findall("ID=(\w+),", line)
        if len(tag) > 0:
          tagType = ''
          if line.count("FORMAT") > 0:
            tagType = 'format'
          elif line.count("INFO") > 0:
            tagType = 'info'
          elif line.count("FILTER") > 0:
            result['filters'][re.findall("ID=(\w+),")[0]] = re.findall('Description="(.*)"')[0]

          typ = re.findall("Type=(\w+),", line)
          if tagType != '':
            number = re.findall("Number=(\w+)", line)
            description = re.findall('Description="(.*)"', line)
            if len(number) == 0:
              number = ['.']
            if len(description) == 0:
              description = ['']
            if "format_"+tag[0] not in elevatedTags:
                result['tags'][tagType][tag[0]] = {'type':typ[0], 'description' : description[0], 'number' : number[0]}
        if line[0] == "#" and line[1] != "#":
          result['columns'] = line.strip()
        if line == '' or line[0] != "#":
            break
    return result

def checkSamContainsRead(samFileName):
    for line in open(samFileName, 'r'):
        if line[0] != "@":
            return True
    return False

def translateTagTypeToColumnType(tag):
  if tag['type'] == "Flag":
    return "boolean"
  if tag['number'] != '1':
    return 'string'
  if tag['type'] == "Integer":
    return 'int32'
  if tag['type'] == "Float":
    return "double"
  return "string"


def mapGatk():

    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:opt/jar/CreateSequenceDictionary.jar'
    print os.environ

    regionFile = open("regions.txt", 'w')
    regionFile.write(job['input']['interval'])

    regionFile.close()

    dxpy.DXFile(job['input']['mappings_file']).wait_on_close()

    dxpy.download_dxfile(job['input']['mappings_file'], "input.bam")
    dxpy.download_dxfile(job['input']['reference_file'], "ref.fa")

    gatkIntervals = open("regions.interval_list", 'w')
    for x in re.findall("-L ([^:]*):(\d+)-(\d+)", job['input']['interval']):
        gatkIntervals.write(x[0] + ":" + x[1] + "-" + x[2] + "\n")
    gatkIntervals.close()

    print "Converting Table to SAM"
    if 1:
        print "Indexing"
        subprocess.check_call("samtools index input.bam", shell=True)
        print "Indexing Reference"
        subprocess.check_call("samtools faidx ref.fa", shell=True)
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.CreateSequenceDictionary REFERENCE=ref.fa OUTPUT=ref.dict" ,shell=True)

        command = job['input']['command'] + job['input']['interval']
        #print command
        print "In GATK"
        subprocess.check_call(command, shell=True)
        #command += " | "

        command = "dx_vcfToVariants2 --table_id %s --vcf_file output.vcf --region_file regions.txt" % (job['input']['tableId'])
        if job['input']['compress_reference']:
            command += " --compress_reference"
        if job['input']['infer_no_call']:
            command += " --infer_no_call"
        if job['input']['compress_no_call']:
            command += " --compress_no_call"

        print "Parsing Variants"
        subprocess.check_call(command, shell=True)

    #job['output']['id'] = job['input']['mappings_table_id']
    job['output']['ok'] = True
