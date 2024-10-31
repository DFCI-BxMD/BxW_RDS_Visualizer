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

if [ "$output_folder" == "/" ];
then
   vm_output="/home/dnanexus/$projname"
   mkdir -p "$vm_output"
else
   parent_dir=$(dirname "$output_folder")
   vm_output="/home/dnanexus/$projname$output_folder"
   mkdir -p "/home/dnanexus/$projname$parent_dir"
   dx download $projName:$output_folder/ -o "/home/dnanexus/$projname$parent_dir" -r -f
fi


rds=${rds_file:-}
if [ -z ${rds} ];
then 
   echo "No data image provided"
else
   fileid=`echo ${rds} | awk -F":" '{print $2}' | sed 's/}//g' | tr -d ' ' | tr -d '"'`
   mkdir -p "${vm_output}/rds_vizualizer" 
   filename=`dx describe $fileid | grep "Name" | awk '{sub(/[^ ]+[ ]+/,"")}1' | sed 's, ,\\ ,g'`
   dx download $fileid -o "${vm_output}/rds_vizualizer/${filename}" -f
   rds_fullpath="${vm_output}/rds_vizualizer/${filename}"
   chmod +777 "$rds_fullpath"
fi

chmod -R +777 "$vm_output"

# start RDS Visualizer ...
 cd "${vm_output}/rds_vizualizer"
 url=https://raw.githubusercontent.com/DFCI-BxMD/BxW_RDS_Visualizer/refs/heads/main/Rizky%20Version/Webapp_source
 wget -P "${vm_output}/rds_vizualizer" $url/DESCRIPTION $url/server.R $url/ui.R
 
 # pull and run Shiny Server docker image
 dx download $DX_PROJECT_CONTEXT_ID:/rk_shiny/rds_vis_maria.tar.gz
 docker load -i rds_vis_maria.tar.gz
 
 # attach our rds_vis app's folder as a volume
 docker run --rm -p 443:3838 -e ROOT=TRUE -e RUNROOTLESS=FALSE -e DISABLE_AUTH=true -e "WORKING_DIR=$vm_output" -v $PWD/rds_visualizer:/srv/shiny-server/ -v /var/run/docker.sock:/var/run/docker.sock -v /home/dnanexus:/home/dnanexus -w /home/dnanexus rds_vis_maria

}

