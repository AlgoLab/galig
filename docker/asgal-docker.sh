#!/bin/bash

# Add local user
# with the same owner as /data
USER_ID=$(stat -c %u /data)
GROUP_ID=$(stat -c %g /data)

echo "Starting with UID:GID $USER_ID:$GROUP_ID"
# We create user:group with the correct uid:gid
groupadd -g "$GROUP_ID" group
useradd --shell /bin/bash -u "$USER_ID" -g group -o -c "" -m user
export HOME=/
chown --recursive "$USER_ID":"$GROUP_ID" /data

if [ -s /data/transcripts.fa ]
then
    if [ -s /data/sample_2.fa ]
    then
    gosu user:group \
	/galig/asgal --multi \
		     -g /data/genome.fa \
		     -a /data/annotation.gtf \
		     -s /data/sample_1.fa \
		     -s2 /data/sample_2.fa \
		     -t /data/transcripts.fa \
		     -o /data/output
    else
    gosu user:group \
	/galig/asgal --multi \
		     -g /data/genome.fa \
		     -a /data/annotation.gtf \
		     -s /data/sample_1.fa \
		     -t /data/transcripts.fa \
		     -o /data/output
    fi
else
    gosu user:group \
    /galig/asgal -g /data/genome.fa \
		 -a /data/annotation.gtf \
		 -s /data/sample_1.fa \
		 -o /data/output
fi
