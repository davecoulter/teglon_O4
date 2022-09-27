CREATE USER 'teglon'@'%' IDENTIFIED WITH mysql_native_password BY '4tegl0n123!!';
GRANT ALL PRIVILEGES ON teglon.* TO 'teglon'@'%' WITH GRANT OPTION;

CREATE USER 'teglon'@'localhost' IDENTIFIED WITH mysql_native_password BY '4tegl0n123!!';
GRANT ALL PRIVILEGES ON YSE.* TO 'teglon'@'localhost' WITH GRANT OPTION;

FLUSH PRIVILEGES;