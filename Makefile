docker: down
	docker build src/ -t gsdb
	docker-compose up -d

down:
	docker-compose down