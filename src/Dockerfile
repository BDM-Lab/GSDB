FROM hekademeia/cheng-common:latest

COPY ./nginx.conf /etc/nginx/http.d/default.conf
COPY ./GSDB /var/www/
RUN echo 'upload_max_filesize = 40M' >> /etc/phpcat 7/php.ini \
  && echo 'post_max_size = 40M' >> /etc/php7/php.ini