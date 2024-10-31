#!/bin/bash

main() {
set -euxo pipefail

# mount the project via dxFuse
mountpoint=$HOME/project
projName=$DX_PROJECT_CONTEXT_ID
wget https://github.com/dnanexus/dxfuse/releases/download/v1.2.0/dxfuse-linux
chmod +x dxfuse-linux
source environment >& /dev/null
echo "Mounting dxfuse"
mkdir -p "$mountpoint"
sudo -E ./dxfuse-linux -uid $(id -u) -gid $(id -g) -verbose 2 -limitedWrite "$mountpoint" "$projName"
projname=`dx describe $DX_PROJECT_CONTEXT_ID | grep "Name" | awk '{sub(/[^ ]+[ ]+/,"")}1' | sed 's, ,\\ ,g'`
output_folder=`dx describe $DX_JOB_ID --json | jq .folder | tr -d '"'`

vm_output="/home/dnanexus"
rstudio_sync_dir=""
if [ "$output_folder" == "/" ];
then
   rstudio_folder=rstudio_$(date +%Y%m%d_%H%M%S)
   vm_output="/home/dnanexus/$projname/$rstudio_folder"
   rstudio_sync_dir="/home/dnanexus/project/$projname/$rstudio_folder"
   mkdir -p "$rstudio_sync_dir"
   mkdir -p "$vm_output"
else
   parent_dir=$(dirname "$output_folder")
   vm_output="/home/dnanexus/$projname$output_folder"
   rstudio_sync_dir="/home/dnanexus/project/$projname$output_folder"
   mkdir -p "/home/dnanexus/$projname$parent_dir"
   dx download $projName:$output_folder/ -o "/home/dnanexus/$projname$parent_dir" -r -f
fi


rdata=${rdata_image:-}
if [ -z ${rdata} ];
then 
   echo "No data image provided"
else
   fileid=`echo ${rdata} | awk -F":" '{print $2}' | sed 's/}//g' | tr -d ' ' | tr -d '"'`
   mkdir -p "${vm_output}/rdata_image" 
   filename=`dx describe $fileid | grep "Name" | awk '{sub(/[^ ]+[ ]+/,"")}1' | sed 's, ,\\ ,g'`
   dx download $fileid -o "${vm_output}/rdata_image/${filename}" -f
   rdata_fullpath="${vm_output}/rdata_image/${filename}"
   chmod +777 "$rdata_fullpath"
fi

chmod -R +777 "$vm_output"

# start RStudio Server ...
docker=${docker_image:-}
if [ -z "${docker}" ];
then
   docker run --rm -p 443:8787 -e ROOT=TRUE -e RUNROOTLESS=FALSE -e DISABLE_AUTH=true -e "WORKING_DIR=$vm_output" -e "PROJECT_DIR=$rstudio_sync_dir" -v /var/run/docker.sock:/var/run/docker.sock -v /home/dnanexus:/home/dnanexus -w /home/dnanexus tariship/rstudio_docker_4.4.0:latest
else
   dockerid=`echo ${docker} | awk -F":" '{print $2}' | sed 's/}//g' | tr -d ' ' | tr -d '"'`
   dx download $dockerid -o /home/dnanexus 
   docker_name_ext=`dx describe $dockerid | grep "Name" | awk '{sub(/[^ ]+[ ]+/,"")}1' | sed 's, ,\\ ,g'`
   docker load < "/home/dnanexus/${docker_name_ext}"
   docker_name=${docker_name_ext%%.tar*}
   docker_name=${docker_name//[-]/:}
   docker run --rm -p 443:8787 -e ROOT=TRUE -e RUNROOTLESS=FALSE -e DISABLE_AUTH=true -e "WORKING_DIR=$vm_output" -e "PROJECT_DIR=$rstudio_sync_dir" -v /var/run/docker.sock:/var/run/docker.sock -v /home/dnanexus:/home/dnanexus -w /home/dnanexus $docker_name
fi
}

