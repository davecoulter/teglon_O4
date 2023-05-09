from celery import Celery
from celery.utils.log import get_task_logger

from web.src.services.teglon import *

logger = get_task_logger(__name__)

app = Celery('teglon_worker_tasks',
             broker='amqp://admin:mypass@rabbit:5672',
             backend='rpc://')

@app.task()
def load_map_batch_process(superevent_id):
    logger.info('Got Request - Starting work ')

    teglon = Teglon()

    # # First load a map...
    print("loading map `%s`..." % superevent_id)
    teglon.load_map(gw_id=superevent_id, clobber=True)
    print("... Done loading `%s`!!" % superevent_id)

    print("extracting map `%s`..." % superevent_id)
    teglon.extract_tiles(gw_id=superevent_id)
    print("... Done extracting `%s`!!" % superevent_id)

    print("plotting map `%s`..." % superevent_id)
    teglon.plot_teglon(gw_id=superevent_id)
    print("... Done plotting `%s`!!" % superevent_id)

    logger.info('Work Finished ')
    return 0