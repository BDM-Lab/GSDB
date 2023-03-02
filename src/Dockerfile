FROM hekademeia/cheng-common:latest

COPY ./nginx.conf /etc/nginx/http.d/default.conf
ADD ./GSDB /var/www/
RUN echo 'upload_max_filesize = 40M' >> /etc/php7/php.ini \
  && echo 'post_max_size = 40M' >> /etc/php7/php.ini