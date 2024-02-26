from flask import Flask, render_template, request, send_from_directory, jsonify, url_for, send_file, Blueprint
import os, io
#import fDeepDNAshape # Assuming you have a Python module named deepDNAshape
import json
import Pyro4
import logging

#logging.basicConfig(level=logging.DEBUG)
#Pyro4.config.DETAILED_TRACEBACK = True

app = Flask(__name__)
#my_app_blueprint = Blueprint('deepdnashape', __name__, url_prefix='/deepdnashape', template_folder = "/srv/www/deepdnashape/templates")

#fDeepDNAshape_predictor = fDeepDNAshape.predictor()
fDeepDNAshape_predictor = Pyro4.Proxy("PYRONAME:deepdnashape.db")
UPLOAD_FOLDER = "/srv/www/deepdnashape/uploads"
DOWNLOAD_FOLDER = "/srv/www/deepdnashape/downloads"
#if not os.path.exists(UPLOAD_FOLDER):
#   os.makedirs(UPLOAD_FOLDER)
#if not os.path.exists(DOWNLOAD_FOLDER):
#    os.makedirs(DOWNLOAD_FOLDER)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['DOWNLOAD_FOLDER'] = DOWNLOAD_FOLDER

nav_items = [
    {"name": "Home", "endpoint": "home"},
    {"name": "Manual", "endpoint": "manual"},
    # Add more items as needed
]

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html', nav_items=nav_items)

@app.route('/')
def home():
    return render_template('index.html', nav_items=nav_items)


# Load manual content from JSON file
with open('/srv/www/deepdnashape/manual/manual.json', 'r') as file:
    manual_content = json.load(file)
@app.route('/manual')
def manual():
    return render_template('manual.html', nav_items=nav_items, content=manual_content)

# Show images
@app.route('/images/<filename>')
def images(filename):
    # Define the directory where your images are stored
    image_directory = 'images'
    # Send the requested image file from the specified directory
    return send_from_directory(image_directory, filename)


interbasepair_features = {"Shift", "Slide", "Rise", "Tilt", "Roll", "HelT", "Shift-FL", "Slide-FL", "Rise-FL", "Tilt-FL", "Roll-FL", "HelT-FL"}
feature_unit = {"Shift": " (Å)", "Slide": " (Å)", "Rise": " (Å)", "Shear": " (Å)", "Stretch": " (Å)", "Stagger": " (Å)", "MGW": " (Å)",
                "Tilt": " (°)", "Roll": " (°)", "HelT": " (°)", "Buckle": " (°)", "ProT": " (°)", "Opening": " (°)",
                "EP": " (kT/e)"}
@app.route('/get_data', methods=['POST'])
def get_data():
    sequences = request.form['sequence'].splitlines()
    selector_value = int(request.form['selector'])
    shape_feature = request.form['shape-features']
    fl_selector = request.form.get("shape-fluctuations")
    k = 2 if shape_feature in interbasepair_features else 1
    #if fl_selector:
    #    shape_feature = shape_feature + "-FL"
    # Add any other processing or arguments as needed
    x_s = []
    y_s = []
    error_ys = []
    legends = []
    texts = []
    for i in range(len(sequences)):
        x = sequences[i]
        y = fDeepDNAshape_predictor.predictSeq(seq = sequences[i], feature = shape_feature, layer = selector_value)
        if fl_selector:
            error_y = fDeepDNAshape_predictor.predictSeq(seq = sequences[i], feature = shape_feature + "-FL", layer = selector_value)
            error_ys.append(error_y)
        x_s.append(list(range(len(x))))
        y_s.append(y)
        legends.append(x)
        texts.append([x[i:i+k] for i in range(len(x))])
    output = {
        "x": x_s,
        "y": y_s,
        "legends": legends,
        "text": texts,
        "shape_feature": shape_feature + feature_unit[shape_feature],
        "error_y": error_ys
    }
    #print(output)
    return jsonify(output)


@app.route('/upload_file', methods=['POST'])
def upload_file():
    if 'file' in request.files:
        file = request.files['file']
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
        file.save(filepath)
        shape_feature = request.form['shape-features']
        fl_selector = request.form.get("shape-fluctuations")
        if fl_selector:
            shape_feature = shape_feature + "-FL"
        selector_value = int(request.form['selector'])
        # Assuming deepDNAshape.process_file returns a path to the generated predictions file
        predictions_path = fDeepDNAshape_predictor.predictFile(file = filepath, feature = shape_feature, layer = selector_value)
        if predictions_path:
        # Return the download link
            return jsonify({
                'success': True,
                'download_link': url_for('download', filename=predictions_path.split('/')[-1], upload = ".".join(file.filename.split(".")[:-1]), feature = shape_feature)
            })
    return jsonify({'success': False})

@app.route('/downloads/<filename>/<upload>/<feature>', methods=['GET'])
def download(filename, upload, feature):
    file_path = os.path.join(app.config["DOWNLOAD_FOLDER"], filename)

    return_data = io.BytesIO()
    with open(file_path, 'rb') as fo:
        return_data.write(fo.read())
    # (after writing, cursor will be at last byte, so move it to start)
    return_data.seek(0)

    os.remove(file_path)
    #shape_feature = request.args.get('shape-features')
    return send_file(return_data, mimetype='application/txt',
                     download_name = upload + "_" + feature + ".txt", as_attachment = True)
    #return send_from_directory(app.config["DOWNLOAD_FOLDER"], filename)

#app.register_blueprint(my_app_blueprint)
if __name__ == "__main__":
    app.run(debug=True)
