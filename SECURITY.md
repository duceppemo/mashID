# Security

mashID runs local command-line tools on local files. It downloads databases from Figshare over HTTPS
and verifies their MD5 against values stored in the package. It never sends data anywhere.

If you find a security problem, for example an unsafe handling of file paths or a way to make the
downloader write outside its directory, please email duceppemo@gmail.com rather than opening a public
issue. You will get an answer within a week.

Only the latest release is supported with fixes.
