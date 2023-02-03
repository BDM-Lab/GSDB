# Readme

## How to update database

The setup process has been moved to install/ and can now be executed by a docker-container using the following setup.

`docker-compose.yml`

```yml
gsdb_init: # initializes the database
  image: gsdb
  depends_on:
    - mysql
  command: php /var/www/install/initdb.php
  links:
    - mysql:mysql
  env_file:
    - env
```