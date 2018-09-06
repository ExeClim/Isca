#!/bin/bash

# Add local user
# Either use the LOCAL_USER_ID if passed in at runtime or
# fallback

if [ -n "$LOCAL_USER_ID" ]; then
	USER_ID=${LOCAL_USER_ID:-9001}
	GROUP_ID=${LOCAL_GROUP_ID:-9003}

	echo "Starting with UID : $USER_ID:$GROUP_ID"
	groupadd -g $GROUP_ID group
	useradd --shell /bin/bash -u $USER_ID  -g $GROUP_ID -o -c "" -m user
	export HOME=/home/user

	exec /usr/local/bin/gosu user "$@"
else
	# run as the isca user
	exec /usr/local/bin/gosu isca "$@"
fi
