#!/bin/bash

# Add local user
# with the same owner as /data
USER_ID=$(/usr/bin/stat -c %u /data)
GROUP_ID=$(/usr/bin/stat -c %g /data)

# When running as root, use a dummy command to change user
SU_CMD="setpriv --reuid=root --regid=root"

# If the owner is not root, we create user:group with the correct uid:gid
if [ "$USER_ID" != "0" ]
then
    echo "Changing user"
    groupadd -g "$GROUP_ID" group
    useradd --shell /bin/bash -u "$USER_ID" -g "$GROUP_ID" -o -c "" -m user
    chown --recursive "$USER_ID":"$GROUP_ID" /data
    SU_CMD="setpriv --reuid="$USER_ID" --regid="$GROUP_ID" --init-groups"
fi
export HOME=/

echo "Starting with UID:GID $USER_ID:$GROUP_ID"
cd /data

if [ -s /data/transcripts.fa ]
then
    if [ -s /data/sample_1.fq ]
    then
	if [ -s /data/sample_2.fq ]
	then
	    $SU_CMD       \
		/galig/asgal --multi \
		             -g /data/genome.fa \
			     -a /data/annotation.gtf \
			     -s /data/sample_1.fq \
			     -s2 /data/sample_2.fq \
			     -t /data/transcripts.fa \
			     -o /data/output
	else
	    $SU_CMD       \
		/galig/asgal --multi \
		             -g /data/genome.fa \
			     -a /data/annotation.gtf \
			     -s /data/sample_1.fq \
			     -t /data/transcripts.fa \
			     -o /data/output
	fi
    else
	if [ -s /data/sample_2.fa ]
	then
	    $SU_CMD       \
		/galig/asgal --multi \
		             -g /data/genome.fa \
			     -a /data/annotation.gtf \
			     -s /data/sample_1.fa \
			     -s2 /data/sample_2.fa \
			     -t /data/transcripts.fa \
			     -o /data/output
	else
	    $SU_CMD       \
		/galig/asgal --multi \
		             -g /data/genome.fa \
			     -a /data/annotation.gtf \
			     -s /data/sample_1.fa \
			     -t /data/transcripts.fa \
			     -o /data/output
	fi
    fi
else
    if [ -s /data/sample_1.fq ]
    then
	$SU_CMD       \
	    /galig/asgal -g /data/genome.fa \
	                 -a /data/annotation.gtf \
			 -s /data/sample_1.fq \
			 -o /data/output
    else
	$SU_CMD       \
	    /galig/asgal -g /data/genome.fa \
	                 -a /data/annotation.gtf \
			 -s /data/sample_1.fa \
			 -o /data/output
    fi
fi
