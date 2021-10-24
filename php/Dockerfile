# Dockerfile
FROM php:8.0-apache

WORKDIR /app
RUN chmod -R 777 /app

# Parametre d apache
COPY conf/vhost.conf /etc/apache2/sites-available/000-default.conf
COPY conf/apache.conf /etc/apache2/conf-available/z-app.conf
COPY conf/php.ini /usr/local/etc/php/conf.d/app.ini
RUN a2enconf z-app

COPY index.php /app/index.php