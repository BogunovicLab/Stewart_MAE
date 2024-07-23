#!/bin/bash

## cp /sc/arion/projects/ISDS/ostewart/MAE_merge/fastqDir/fastq_filt/auto_trim.sh /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/scripts
## cp /sc/arion/projects/ISDS/ostewart/MAE_merge/fastqDir/fastq_filt/auto_trim.lsf /sc/arion/projects/ISDS/haleyr/RMAE_OJAY/scripts
## bsub < auto_trim.lsf

# Usage
# This script automates running fastp for multiple PE data.
# Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>

# Execution
# Case(1) run on a couple of PE files with extension *.fq
# $ source auto_trim.sh *.fq

# Case(2) run on a couple of PE files with extension *.fq.gz
# $ source auto_trim.sh *.fq.gz

# Case(3) run on a all the fastq files in the current directory (mixed extensions, like .fq. .fastq )
# $ source auto_trim.sh *

# Invoke help
# $ source auto_trim.sh -h
# $ source auto_trim.sh --help
# Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>

ml fastp

red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
reset=`tput sgr0`
count=0

usage ()
{
  echo -e "${green}Usage: source auto_trim.sh [*.extension]\n \
      extension: <fq> or <fastq> or <fq.gz> or <fastq.gz>\n \
      example: source auto_trim.sh *.fq.gz\n ${reset}\n\
${yellow}Help:  source auto_trim.sh -h or --help${reset}"
  return
}

file_not_found ()
{
echo -e "\n${red}FileNotFoundError: No such file with extension $@ found!${reset}"
echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}\n"
return 
}

file_name_error ()
{
echo -e "\n${red}Filename Error: Paired end file names should contain _R1 _R2${reset}"
echo -e "${green}Example: test_R1.fq.gz, test_R2.fq.gz${reset}\n"
return 
}

file_extension_error ()
{
echo -e "\n${red}FileExtensionError: Invalid extension${reset}"
echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}\n"
return     
}

if [[ ( $1 == '-h' ) || ( $1 == '--help') ]] ;then
 usage
elif [[ $# -eq 0 ]] ;then
 echo "${red}Error: No parameter(s) provided${reset}"
 usage
 return 0
else
  for i in $@; do
        count=$((count+1))
        if [ -f $i ] ;then
            if [[ (${i#*.} == "fastq.gz") || (${i#*.} == "fq.gz") || (${i#*.} == "fastq") || (${i#*.} == "fq") ]] ;then
            	 if echo $1 | grep -q -e "_R1" -e "_R2"; then
                   if [[ $count%2 -ne 0 ]]; then
                       sample_name=`echo $i | awk -F "_R1"  '{print $1}'`
                       extension=`echo $i | awk -F "_R1"  '{print $2}'`
                       R1=${sample_name}_R1${extension}
                       R2=${sample_name}_R2${extension}
                       R1_filt=${sample_name}_R1_filt${extension}
                       R2_filt=${sample_name}_R2_filt${extension}
                       echo -e "\n${yellow}[Running fastp for sample] ${sample_name} at `whoami`${reset}\n"
                       date && time fastp -i $R1 -I $R2 \
                        --verbose --detect_adapter_for_pe --overrepresentation_analysis \
                        --thread 8 \
                        -o $R1_filt -O $R2_filt

                   fi
                else file_name_error
                fi      
            elif [[ (${i#*.} == "sh") || (${i#*.} == "sh~")  ]] ;then
                 echo -n
            else
                echo -ne "${red}Check:$i${reset}"
                file_extension_error  
            fi     
        else
            file_not_found $@
        fi
  done
fi