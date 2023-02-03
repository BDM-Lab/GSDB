docker: down
	docker build src/ -t gsdb
	docker-compose up -d

down:
	docker-compose down

publish:
#	docker build src/ -t hekademeia/gsdb:latest