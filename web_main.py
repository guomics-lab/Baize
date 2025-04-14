from flask import Flask

from src.baize.baize_router import baize_bp



app = Flask(__name__)

class CaseInsensitiveMiddleware:
    def __init__(self, app):
        self.app = app

    def __call__(self, environ, start_response):
        path_info = environ['PATH_INFO']
        lower_path = path_info.lower()
        environ['PATH_INFO'] = lower_path

        return self.app(environ, start_response)

app.wsgi_app = CaseInsensitiveMiddleware(app.wsgi_app)

app.register_blueprint(baize_bp, url_prefix='/baize')


if __name__ == '__main__':
    app.run(debug=True, port=9997, host='0.0.0.0')

