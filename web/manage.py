from flask.cli import FlaskGroup
import os

# print("hello world")
# print(os.path.dirname(__file__))
from web.project import app
# from project import app


cli = FlaskGroup(app)


if __name__ == "__main__":
    cli()