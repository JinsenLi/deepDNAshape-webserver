import sys
import logging
 
sys.path.insert(0, '/srv/www/deepdnashape')
sys.path.insert(0, '/home/jinsenli/deepdnashape-env/lib/python3.9/site-packages/')
 
# Set up logging
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
 
# Import and run the Flask app
from app import app as application