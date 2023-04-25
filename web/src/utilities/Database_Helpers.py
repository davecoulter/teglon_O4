import mysql.connector
from mysql.connector import Error
import MySQLdb as my
# import mysqlclient as my
from configparser import RawConfigParser
import time
import sys
# import os
# print(os.path.dirname(__file__))

# configFile = '/opt/project/Settings.ini'
configFile = 'Settings.ini'
# configFile = '../Settings.ini'
config = RawConfigParser()
config.read(configFile)

db_name = config.get('database', 'DATABASE_NAME')
db_user = config.get('database', 'DATABASE_USER')
db_pwd = config.get('database', 'DATABASE_PASSWORD')
db_host = config.get('database', 'DATABASE_HOST')
db_port = config.get('database', 'DATABASE_PORT')

def bulk_upload(query):
    success = False
    try:
        conn = mysql.connector.connect(user=db_user, password=db_pwd, host=db_host, port=db_port, database=db_name,
                                       allow_local_infile=True)
        cursor = conn.cursor()
        cursor.execute(query)
        conn.commit()
        success = True

    except Error as e:
        print("Error in uploading CSV!")
        print(e)
    except Exception as e:
        print("Error in uploading CSV!")
        print(e)
        test = 1
    finally:
        cursor.close()
        conn.close()

    return success

def query_db(query_list, commit=False):

    results = []
    try:
        chunk_size = 1e+6
        db = my.connect(host=db_host, user=db_user, passwd=db_pwd, db=db_name, port=3306)
        cursor = db.cursor()

        query_count = len(query_list)
        for qi, q in enumerate(query_list):
            current_q = qi + 1
            print("Executing %s/%s" % (current_q, query_count))
            cursor.execute(q)

            if commit:  # used for updates, etc
                db.commit()

            streamed_results = []
            print("\tfetching results...")
            while True:
                r = cursor.fetchmany(1000000)
                count = len(r)
                streamed_results += r
                size_in_mb = sys.getsizeof(streamed_results) / 1.0e+6

                print("\t\tfetched: %s; current length: %s; running size: %0.3f MB" % (
                count, len(streamed_results), size_in_mb))

                if not r or count < chunk_size:
                    break

            results.append(streamed_results)

    # except Error as e:
    except Exception as e:
        print('Exception:', e)
    finally:
        cursor.close()
        db.close()

    print("Returning total results: %s" % len(results))
    return results

def batch_query(query_list):
    return_data = []
    batch_size = 500
    ii = 0
    jj = batch_size
    kk = len(query_list)

    print("\nLength of data to query: %s" % kk)
    print("Query batch size: %s" % batch_size)
    print("Starting loop...")

    number_of_queries = len(query_list) // batch_size
    if len(query_list) % batch_size > 0:
        number_of_queries += 1

    query_num = 1
    payload = []
    while jj < kk:
        t1 = time.time()

        print("%s:%s" % (ii, jj))
        payload = query_list[ii:jj]
        return_data += query_db(payload)

        ii = jj
        jj += batch_size
        t2 = time.time()

        print("\n********* start DEBUG ***********")
        print("Query %s/%s complete - execution time: %s" % (query_num, number_of_queries, (t2 - t1)))
        print("********* end DEBUG ***********\n")

        query_num += 1

    print("Out of loop...")

    t1 = time.time()

    print("\n%s:%s" % (ii, kk))

    payload = query_list[ii:kk]
    return_data += query_db(payload)

    t2 = time.time()

    print("\n********* start DEBUG ***********")
    print("Query %s/%s complete - execution time: %s" % (query_num, number_of_queries, (t2 - t1)))
    print("********* end DEBUG ***********\n")

    return return_data

def insert_records(query, data):
    _tstart = time.time()
    success = False
    try:
        conn = mysql.connector.connect(user=db_user, password=db_pwd, host=db_host, port=db_port, database=db_name,
                                       allow_local_infile=True)
        cursor = conn.cursor()
        cursor.executemany(query, data)

        conn.commit()
        success = True
    except Error as e:
        print('Error:', e)
    finally:
        cursor.close()
        conn.close()

    _tend = time.time()
    print("\n********* start DEBUG ***********")
    print("insert_records execution time: %s" % (_tend - _tstart))
    print("********* end DEBUG ***********\n")
    return success

def batch_insert(insert_statement, insert_data, batch_size=50000):
    _tstart = time.time()

    i = 0
    j = batch_size
    k = len(insert_data)

    print("\nLength of data to insert: %s" % len(insert_data))
    print("Insert batch size: %s" % batch_size)
    print("Starting loop...")

    number_of_inserts = len(insert_data) // batch_size
    if len(insert_data) % batch_size > 0:
        number_of_inserts += 1

    insert_num = 1
    payload = []
    while j < k:
        t1 = time.time()

        print("%s:%s" % (i, j))
        payload = insert_data[i:j]

        if insert_records(insert_statement, payload):
            i = j
            j += batch_size
        else:
            raise ("Error inserting batch! Exiting...")

        t2 = time.time()

        print("\n********* start DEBUG ***********")
        print("INSERT %s/%s complete - execution time: %s" % (insert_num, number_of_inserts, (t2 - t1)))
        print("********* end DEBUG ***********\n")

        insert_num += 1

    print("Out of loop...")

    t1 = time.time()

    print("\n%s:%s" % (i, k))

    payload = insert_data[i:k]
    if not insert_records(insert_statement, payload):
        raise ("Error inserting batch! Exiting...")

    t2 = time.time()

    print("\n********* start DEBUG ***********")
    print("INSERT %s/%s complete - execution time: %s" % (insert_num, number_of_inserts, (t2 - t1)))
    print("********* end DEBUG ***********\n")

    _tend = time.time()

    print("\n********* start DEBUG ***********")
    print("batch_insert execution time: %s" % (_tend - _tstart))
    print("********* end DEBUG ***********\n")

def delete_rows(delete_query, id_tuples_to_delete):
    # Start with
    batch_size = 10000
    i = 0
    j = batch_size
    k = len(id_tuples_to_delete)

    print("\nLength of records to DELETE: %s" % k)
    print("DELETE batch size: %s" % j)
    print("Starting loop...")

    number_of_deletes = k // batch_size
    if k % batch_size > 0:
        number_of_deletes += 1

    delete_num = 1
    while j < k:
        t1 = time.time()
        print("%s:%s" % (i, j))

        id_string = ",".join([str(id_tup[0]) for id_tup in id_tuples_to_delete[i:j]])
        query_db([delete_query % id_string], commit=True)
        i = j
        j += batch_size

        t2 = time.time()

        print("\n********* start DEBUG ***********")
        print("DELETE %s/%s complete - execution time: %s" % (delete_num, number_of_deletes, (t2 - t1)))
        print("********* end DEBUG ***********\n")

        delete_num += 1

    t1 = time.time()
    print("%s:%s" % (i, k))

    id_string = ",".join([str(id_tup[0]) for id_tup in id_tuples_to_delete[i:k]])
    query_db([delete_query % id_string], commit=True)

    t2 = time.time()

    print("\n********* start DEBUG ***********")
    print("DELETE %s/%s complete - execution time: %s" % (delete_num, number_of_deletes, (t2 - t1)))
    print("********* end DEBUG ***********\n")
