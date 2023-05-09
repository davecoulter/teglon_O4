from flask import Flask, jsonify, send_from_directory
from flask import request
from celery import Celery

app = Flask(__name__)

teglon_worker = Celery('teglon_worker',
                    broker='amqp://admin:mypass@rabbit:5672',
                    backend='rpc://')


@app.route('/load_map/<superevent_id>')
def call_method(superevent_id):
    app.logger.info("Invoking load_map_batch_process()")
    r = teglon_worker.send_task('web.src.services.teglon_worker_tasks.load_map_batch_process',
                                kwargs={'superevent_id': superevent_id})
    app.logger.info(r.backend)
    return r.id


@app.route('/load_status/<task_id>')
def get_status(task_id):
    status = teglon_worker.AsyncResult(task_id, app=teglon_worker)
    print("Invoking Method ")
    return "Status of the Task " + str(status.state)


# @app.route('/simple_task_result/<task_id>')
# def task_result(task_id):
#     result = simple_app.AsyncResult(task_id).result
#     return "Result of the Task " + str(result)







@app.route("/")
def hello_world():
    test = 1
    return jsonify(hello="world")





# @app.route("/")
# def load_map():
#     test = 1
#
#     # 1. Do work
#     # 2. What to return?
#
#     pass
#     # return jsonify(hello="world")
#
# @app.route("/static/<path:filename>")
# def staticfiles(path, filename):
#     test = 1
#     return send_from_directory("/app/web/events/%s" % path, filename)
#
# @app.route("/get_tiles", methods=["GET"])
# def get_tiles():
#     # def extract_tiles(gw_id, healpix_file, s_band='r', t_band='r', g_band='r', s_exp_time=60.0, t_exp_time=120.0,
#     #                   g_exp_time=120.0, extinct=0.5, prob_type='4D', s_cum_prob=0.9, t_cum_prob=0.9, s_min_ra=-1,
#     #                   s_max_ra=-1, s_min_dec=-1, s_max_dec=-1, t_min_ra=-1, t_max_ra=-1, t_min_dec=-1, t_max_dec=-1,
#     #                   g_min_ra=-1, g_max_ra=-1, g_min_dec=-1, g_max_dec=-1, galaxies_detector='NICKEL', num_gal=10000):
#
#     telescope = request.args.get('telescope')
#
#     gw_id = request.args.get('gw_id')
#     healpix_file = request.args.get('healpix_file')
#     prob_type = request.args.get('prob_type')
#     try:
#         s_cum_prob = float(request.args.get('s_cum_prob'))
#         t_cum_prob = float(request.args.get('t_cum_prob'))
#     except:
#         return jsonify({
#             gw_id: gw_id,
#             healpix_file: healpix_file,
#             telescope: telescope,
#             msg: "Wrong s_cum_prob/t_cum_prob type (must be a float)"
#         })
#
#     csv_output_path = extract_tiles(gw_id=gw_id, healpix_file=healpix_file, prob_type=prob_type, s_cum_prob=s_cum_prob,
#                   t_cum_prob=t_cum_prob)
#
#     output_tokens = csv_output_path.split("/")
#
#     return jsonify({
#         "gw_id": gw_id,
#         "healpix_file": healpix_file,
#         "telescope": telescope,
#         "file_output": "http://0.0.0.0:1337/static/%s/%s" % (output_tokens[-2], output_tokens[-1])
#     })