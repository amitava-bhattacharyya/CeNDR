import csv
import logging
from flask import Flask
from flask_restful import Api
from flask_debugtoolbar import DebugToolbarExtension
from models import *
from datetime import date, datetime
from urlparse import urljoin
from flask.ext.cache import Cache

# Fetch credentials
from gcloud import datastore
ds = datastore.Client(project="andersen-lab")

def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""

    if isinstance(obj, date):
        serial = obj.isoformat()
        return serial
    raise TypeError("Type not serializable")


def autoconvert(s):
    for fn in (int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

# Caching
app = Flask(__name__, static_url_path='/static')

# Cache
cache = Cache(app, config={'CACHE_TYPE': 'gaememcached'})

api = Api(app)
app.config['SECRET_KEY'] = os.getenv("SECRET_KEY")
app.config['DEBUG_TB_INTERCEPT_REDIRECTS'] = False


if os.getenv('SERVER_SOFTWARE') and \
        os.getenv('SERVER_SOFTWARE').startswith('Google App Engine/'):
    app.debug = False
    from flask_sslify import SSLify
    sslify = SSLify(app, skips=['strains/global-strain-map'])
else:
    app.debug = True
    app.config['SECRET_KEY'] = "test"
    toolbar = DebugToolbarExtension(app)


#
# Utility Functions
# 

def get_stripe_keys():
    # Retrieves stripe keys
    if (os.getenv('SERVER_SOFTWARE') and
            os.getenv('SERVER_SOFTWARE').startswith('Google App Engine/')):
        return ds.get(ds.key("credential", "stripe_live"))
    else:
        return ds.get(ds.key("credential", "stripe_test"))


from views import *
from cegwas import *