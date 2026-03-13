## Core Unix commands

Most unix commands have super short names, which makes them quick
to type but annoying to learn. Major file system related
commands include.

pwd: Print working directory
ls: List files and directories
cd: Change Directories
mkdir: Make directories
rm: remove files and directories (delete)
cp: copy files (source > destination)
mv: Move files (basically re-name)
nano: A wee text editor (very basic but always available)

curl:     Download files
wget:     Download files from the web
tar -zxvf:UnTar (unpackage) Tar archive files
gunzip :  UnZip files
$PATH:    The places (dirs) to look for programs

## AWS EC2 Instance

Connect to my instance with:
ssh -i ~/Downloads/bimm143_eloise.pem ubuntu@ec2-54-202-63-98.us-west-2.compute.amazonaws.com

Github link:
https://github.com/EloiseSimpson/Bimm143_githhub/tree/main

Secure copy files between machines, in this case from our instance to our laptop
scp -i ~/Downloads/bimm143_eloise.pem ubuntu@ec2-54-202-63-98.us-west-2.compute.amazonaws.com:/home/ubuntu/work/Bimm143_githhub/Class16/results.txt .

## Class 17 instance
ssh -i "bimm143_eloise.pem" ubuntu@ec2-52-33-145-234.us-west-2.compute.amazonaws.com

export KEY=~/Downloads/bimm143_eloise.pem
export SERVER=ubuntu@ec2-52-33-145-234.us-west-2.compute.amazonaws.com

ssh -i $KEY $SERVER
scp -r 


