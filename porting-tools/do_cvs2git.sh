#! /bin/sh


# in the cvs2git27 conda env
cvs2git \
  --use-rcs \
  --blobfile=/users/livesey/emls/mlspgs-git-migration-17jun2025/cvs2git-out/mlspgs.blob \
  --dumpfile=/users/livesey/emls/mlspgs-git-migration-17jun2025/cvs2git-out/mlspgs.git-dump \
  --tmpdir=/users/livesey/emls/mlspgs-git-migration-17jun2025/cvs2git-tmp \
  --encoding=utf_8 \
  --encoding=iso-8859-1 \
  --encoding=windows-1252 \
  /users/livesey/cvsroot-mlsdev-backup/cvsroot/mlspgs
