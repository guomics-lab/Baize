import hashlib
import logging
import os
import shutil

from flask import render_template, request, send_file, send_from_directory, Blueprint

from src.baize import containmination_index_calculation

baize_bp = Blueprint('Baize', __name__)

base_save_dir = './upload_data'

demo_data_dir = './demo_data'

logger = logging.getLogger()
logger.setLevel(logging.INFO)

logFile = './baize.log'
fh = logging.FileHandler(logFile, mode='a', encoding='utf-8')
fh.setLevel(logging.INFO)
format = logging.Formatter('%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s')
fh.setFormatter(format)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(format)
logger.addHandler(ch)


def get_str_md5(content):
    """
    :param content:
    :return:
    """
    m = hashlib.md5(content.encode())
    return m.hexdigest()


@baize_bp.route('/static/<path:filename>')
def static_file(filename):
    return send_from_directory('static/baize', filename)


@baize_bp.route('', methods=['GET'])
def index():
    return render_template('baize/index.html')


@baize_bp.route('/', methods=['GET'])
def index2():
    return render_template('baize/index.html')


@baize_bp.route('/uploadsample', methods=['POST'])
def upload_sample():
    try:
        if 'file' not in request.files:
            logger.error('No file part')
            return 'Parameter error'
        file = request.files['file']
        pageUUID = request.form.get('pageUUID')
        if not pageUUID:
            logger.error('No pageUUID')
            return 'Parameter error'
        if file.filename == '':
            logger.error('file is empty')
            return 'Parameter error'
        logger.info(f'upload_sample: pageUUID: {pageUUID}, fileName: {file.filename}')
        if file and allowed_file(file.filename):
            filename = file.filename
            save_dir = os.path.join(base_save_dir, get_str_md5(pageUUID))
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)

            file.save(os.path.join(save_dir, filename))
            return 'File successfully uploaded'
        return 'The file type does not meet the requirements'
    except Exception as e:
        logger.exception(f'upload_sample error:')
        return 'Upload failed'


@baize_bp.route('/uploadprotein', methods=['POST'])
def upload_protein():
    try:
        if 'file' not in request.files:
            return 'Parameter error'
        pageUUID = request.form.get('pageUUID')
        if not pageUUID:
            return 'Parameter error'
        file = request.files['file']
        if file.filename == '':
            logger.error('file is empty')
            return 'Parameter error'
        logger.info(f'upload_sample: pageUUID: {pageUUID}, fileName: {file.filename}')
        if file and allowed_file(file.filename):
            save_filename = 'ProteinFile_' + file.filename
            save_dir = os.path.join(base_save_dir, get_str_md5(pageUUID))
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            file.save(os.path.join(save_dir, save_filename))
            return 'File successfully uploaded'
        return 'The file type does not meet the requirements'
    except Exception as e:
        logger.exception(f'upload_protein error:')
        return 'Upload failed'


@baize_bp.route('/calc', methods=['POST'])
def calc():
    try:
        pageUUID = request.form.get('pageUUID')
        proteinFileName = request.form.get('proteinFileName')
        sampleFileName = request.form.get('sampleFileName')
        uniprotIDColumnName = request.form.get('uniprotIDColumnName', '').strip()
        geneColumnName = request.form.get('geneColumnName', '').strip()
        isCoagulationChecked = request.form.get('isCoagulationChecked')
        isRbcChecked = request.form.get('isRbcChecked')
        isPlateletChecked = request.form.get('isPlateletChecked')

        if not pageUUID or not proteinFileName or not sampleFileName:
            return 'Parameter error'

        pageUUID = get_str_md5(pageUUID)

        each_base_save_dir = os.path.join(base_save_dir, pageUUID)
        pg_matrix = os.path.join(each_base_save_dir, 'ProteinFile_' + proteinFileName)
        sampleinfo = os.path.join(each_base_save_dir, sampleFileName)

        contaminationTypeList = []
        if isCoagulationChecked == 'true':
            contaminationTypeList.append('coagulation')
        if isRbcChecked == 'true':
            contaminationTypeList.append('rbc')
        if isPlateletChecked == 'true':
            contaminationTypeList.append('platelet')

        if len(contaminationTypeList) == 0:
            return 'Please select the type of contamination'

        contaminationType = ','.join(contaminationTypeList)

        try:
            containmination_index_calculation.main(each_base_save_dir, pg_matrix, sampleinfo, contaminationType,
                                                   uniprotIDColumnName, geneColumnName)

            return 'File calculation completed'
        except Exception as e:
            logger.exception(f'containmination_index_calculation error:')
            return 'Calculate error'
    except Exception as e:
        logger.exception(f'upload_protein error:')
        return 'Calculate exception'


@baize_bp.route('/img', methods=['GET'])
def display_img():
    pageUUID = request.args.get('pageUUID')
    imgType = request.args.get('imgType')
    pageUUID = get_str_md5(pageUUID)
    each_base_save_dir = os.path.join(base_save_dir, pageUUID, 'result')
    # imgPath = os.path.join(each_base_save_dir, 'OmniProt_contamination_Calculator.png')

    if imgType == 'platelet':
        imgPath = os.path.join(each_base_save_dir, 'Platelet_contamination_index_profile.png')
    elif imgType == 'rbc':
        imgPath = os.path.join(each_base_save_dir, 'Erythrocyte_contamination_index_profile.png')
    elif imgType == 'coagulation':
        imgPath = os.path.join(each_base_save_dir, 'Coagulation_contamination_index_profile.png')

    return send_file(imgPath, mimetype='image/png')


@baize_bp.route('/downfile', methods=['GET'])
def down_file():
    pageUUID = request.args.get('pageUUID')
    pageUUID = get_str_md5(pageUUID)
    project_dir = os.path.join(base_save_dir, pageUUID)
    project_result_dir = os.path.join(project_dir, 'result')
    pp_file_list = os.listdir(project_dir)
    pg_file_name = ''
    for pp_file in pp_file_list:
        if pp_file.startswith('ProteinFile_'):
            pg_file_name = pp_file
            break

    base_pg_file_name = pg_file_name.removeprefix('ProteinFile_')
    base_pg_file_name = os.path.splitext(base_pg_file_name)[0]
    temp_dir_name = f'{base_pg_file_name}_OmniProt_contamination_calculation'
    temp_zip_dir = os.path.join(project_dir, temp_dir_name)
    if os.path.exists(temp_zip_dir):
        shutil.rmtree(temp_zip_dir)
    os.makedirs(temp_zip_dir)
    for file_name in os.listdir(project_result_dir):
        if file_name.endswith('.png'):
            continue
        shutil.copy(os.path.join(project_result_dir, file_name), os.path.join(temp_zip_dir, file_name))

    shutil.make_archive(temp_zip_dir, 'zip', project_dir, temp_dir_name)
    return send_file(os.path.join(project_dir, f'{temp_dir_name}.zip'), as_attachment=True)


@baize_bp.route('/demosample', methods=['GET'])
def demo_sample():
    return send_file(os.path.join(demo_data_dir, 'wosp24266_sampleinfo.xlsx'), as_attachment=True)


@baize_bp.route('/demomatrix', methods=['GET'])
def demo_matrix():
    return send_file(os.path.join(demo_data_dir, 'wosp24266_pg_matrix.csv'), as_attachment=True)


def allowed_file(filename):
    """
    """
    ALLOWED_EXTENSIONS = {'xlsx', 'csv', 'tsv'}
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
