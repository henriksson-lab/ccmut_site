rung:
	gunicorn app:server -b :8052 -n ccmut --timeout=200
install:
	pip install -r requirements.txt
