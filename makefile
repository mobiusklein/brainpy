dev:
	python setup.py develop

test:
	py.test -v  brainpy --cov=brainpy --cov-report=html --cov-report term

retest:
	py.test -v brainpy --lf