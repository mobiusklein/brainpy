dev:
	python setup.py develop

test:
	py.test -v  ./tests --cov=brainpy --cov-report=html --cov-report term

retest:
	py.test -v ./tests --lf

update-docs:
	git checkout gh-pages
	git pull origin master
	cd docs && make html
	git add docs/_build/html -f
	git commit -m "update docs"
	git push origin gh-pages
	git checkout master
