from flask.cli import FlaskGroup

from web.project import app


cli = FlaskGroup(app)


if __name__ == "__main__":
    cli()