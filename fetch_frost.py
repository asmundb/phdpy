import datetime
import requests
import influxdb

def get_frost(start,end,variable,sub_idList=["SN17850"],frost_client_id=""):
    parameters = dict()
    parameters["sources"] = ','.join(sub_idList)
    parameters["elements"] = variable
    parameters['referencetime'] = '%s/%s' % (start, end)
    parameters['timeresolutions'] = "PT1H"
    r = requests.get('https://frost.met.no/observations/v0.jsonld', parameters, auth=(frost_client_id, ''))
    #print(r)
    dat = r.json()["data"]
    #print(dat)
    nlevel = len(dat[0]["observations"])
    nobs = len(dat)
    #sm = ma.zeros((nobs, nlevel))
    levels = [100,5,50,60,75,10,20,30,40]
    levels.sort()
    refTime = []
    data = []
    for it in range(nobs):
        refTime += [dat[it]["referenceTime"]]
        for il in range(nlevel):
            if "level" in dat[it]["observations"][il]:
                lev = dat[it]["observations"][il]["level"]["value"]
            else:
                lev = 0
            swap = False
            if not swap: # manually swap sensors
               sid = dat[it]["sourceId"]
            else:
                sid0 = dat[it]["sourceId"].split(":")
                sid1 = abs(int(sid0[-1])-1)
                sid = "%s:%d" % (sid0[0],sid1)
            data += [{"measurement":"soil",
                      "fields":{variable:float(dat[it]["observations"][il]["value"])},
                      "tags":{"level":"%03d"%int(lev),
                              "sourceId":sid},
                      "time":int(datetime.datetime.strptime(dat[it]["referenceTime"][0:16],"%Y-%m-%dT%H:%M").timestamp())
                      }]
    return data

def update_database(data):
    db = influxdb.InfluxDBClient()
    db.write_points(data,time_precision='s',database='aas')


def reset_db(name):
    db = influxdb.InfluxDBClient()
    db.drop_database(name)
    db.create_database(name)



dt = datetime.datetime.utcnow()
#start = (dt).strftime("%Y-%m-%dT%H")
start = (dt-datetime.timedelta(hours=3)).strftime("%Y-%m-%dT%H")
#start = datetime.datetime(2021,6,9,7).strftime("%Y-%m-%dT%H")
end = dt.strftime("%Y-%m-%dT%H")

vars = ["volume_fraction_of_water_in_soil",
        "soil_temperature",
        "sum(precipitation_amount PT1H)",
        "air_temperature","wind_speed"]#,
#        "mean(relative_humidity)"]
#        "mean(surface_downwelling_shortwave_flux_in_air PT1H)",
#        "mean(surface_downwelling_longwave_flux_in_air PT1H)",
#        "mean(surface_net_downward_radiative_flux PT1H)"]

stlist = ["SN17850:%d" % i for i in (0,1)]

with open("frostId",'r') as f:
  frost_ID = f.readline()



data = []

#var = vars[0]
for var in vars:
    try:
        data += get_frost(start,end,var,sub_idList=stlist, frost_client_id=frost_ID)
    except Exception as e:
        print("failed",var)
        print(e)
if len(data) > 0:
    update_database(data)






