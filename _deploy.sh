#!/bin/sh

#-----------------------------------------------------------------------
# Upload.

cp ./config/logo.png ./_book/config/

echo "Uploading files to server.\n"
rsync -avzp \
      ./_book/ \
      --progress \
      --rsh="ssh -p$PATAXOP" "$PATAXO:~/public_html/cursoR/mrnl/"

#-----------------------------------------------------------------------
