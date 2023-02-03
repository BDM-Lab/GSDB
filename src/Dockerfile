FROM expressiveanalytics/alpine-php7:latest

COPY ./GSDB /var/www/
COPY ./nginx.conf /etc/nginx/http.d/default.conf
RUN apk --no-cache add php7-mysqli php7-pdo php7-pdo_mysql openjdk11