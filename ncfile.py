import netCDF4 as nc
import datetime
import numpy as np
import pyproj
import pyresample
import numpy.ma as ma
import metio.grid
import metio.dataset
import re


class NcFile:
    """
        Read variable, location and time from ncfile
        __call__(self,varname):
            return metio.dataset.Grid object

        get_area_def(self):
            return pyresample.geometry.AreaDefinition object

    """

    name_dict = {"time": "time",
                 "lat": "latitude",
                 "lon": "longitude",
                 "x:": "x",
                 "y:": "y",
                 "proj_vars": "Projection_parameters"
                 }

    def __init__(self, filename):
        self.filename = filename
        if type(filename) == str and "*" not in filename:
            fh = nc.Dataset
        else:
            fh = nc.MFDataset
        self.fh = fh(filename, 'r')
        self.times = _nctime2datetime(self.fh.variables[self.name_dict["time"]])
        self.loc = self._get_loc()

    def __call__(self, varname):
        var = self.fh.variables[varname]
        if hasattr(var, "_FillValue"):
            arr = ma.masked_values(var[:], var._FillValue)
        else:
            arr = var[:]
        return metio.dataset.Grid(self.times, self.loc, arr)

    def close(self):
        self.fh.close()

    def search(self, string):
        out = []
        for v in self.fh.variables:
            if string.lower() in v.lower():
                out.append(v)
        return out

    def get_area_def(self, proj_var="Projection_parameters"):
        pp = self.fh.variables[proj_var]
        try:  # obtain proj4 string
            if hasattr(pp, "proj4"):
                proj4 = pp.proj4
            else:
                if pp.grid_mapping_name == "lambert_conformal_conic":
                    proj_id = "lcc"
                    lat0 = pp.latitude_of_projection_origin
                    lon0 = pp.longitude_of_central_meridian
                    sp = pp.standard_parallel
                    if type(sp) == list and len(sp) > 1:
                        lat1 = sp[0]
                        lat2 = sp[1]
                    else:
                        lat1 = sp
                        lat2 = sp
                    ptmp = "+proj=lcc +lat_1=%.2f +lat_2=%.2f +lat_0=%.2f +lon_0=%.2f +units=m +R=%9.3e +no_defs"
                    proj4 = ptmp % (lat0, lat1, lat2, lon0, pp.earth_radius)
                else:
                    raise NotImplementedError
            proj_id = re.search('\+proj=([a-z]*) \+', proj4).group(1)
            p = pyproj.Proj(proj4, preserve_units=False)
            lons = self.loc.lons
            lats = self.loc.lats
            ny, nx = lons.shape[-2:]
            llx, lly = p(lons[0, 0], lats[0, 0])
            urx, ury = p(lons[-1, -1], lats[-1, -1])
            extent = (llx, lly, urx, ury)
            area_id = "doImatter"
            description = "domain_area"
            return pyresample.geometry.AreaDefinition(area_id, description, proj_id, proj4, nx, ny, extent)
        except Exception as e:
            raise e

    def _get_loc(self):
        lons = self.fh.variables[self.name_dict["lon"]]
        lats = self.fh.variables[self.name_dict["lat"]]
        if len(lons.shape) == 3:
            lons = lons[0, :, :]
            lats = lats[0, :, :]
        loc = metio.grid.Grid(lats, lons)
        return loc


class SfxFa(NcFile):

    pass


class MetNc(NcFile):

    def __init__(self, filename):
        self.name_dict = super().name_dict
        self.name_dict["proj_vars"] = "projection_lcc"
        super().__init__(filename)


def _nctime2datetime(tvar):
    tl = tvar.units.split("since")
    refdate = tl[1].split('+')
    tzone = "0000"
    if len(refdate) > 1:
        tzone = refdate[1].replace(":", "")
    dt = datetime.datetime.strptime(refdate[0].strip() + " +" + tzone, "%Y-%m-%d %H:%M:%S %z")
    if tl[0].strip() == "seconds":
        tout = [dt + datetime.timedelta(seconds=tvar[i].tolist()) for i in range(tvar[:].size)]
    else:
        raise NotImplementedError
    return tout


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

    filename = "/lustre/storeA/users/asmundb/H2O/case_studies/20190713/os3dsd/ref_24h/ICMSHHARM.nc"

    #filename = "/home/asmundb/Projects/H2O/nc/20190713/os3dsd/ref_24h/ICMSHHARM.nc"

    ref = NcFile(filename)
    pp = ref.get_area_def()
    sm = ref("X002WG1")
    mnapath = "/lustre/storeB/project/metkl/klinogrid/archive/met_nordic_analysis/v2/"
    fn2 = [mnapath + "/2019/07/13/met_analysis_1_0km_nordic_20190713T%02dZ.nc" % (i) for i in range(24)]
    #fn2 = ["/home/asmundb/Projects/H2O/nc/nor_ana/met_analysis_1_0km_nordic_20190713T%02dZ.nc" % (i) for i in range(22)]
    ana = MetNc(fn2)
    pp1 = ana.get_area_def(proj_var=ana.search("proj")[0])
    tp = ana("precipitation_amount")
    #tp = ana("air_temperature_2m")

    tp1 = np.moveaxis(pyresample.kd_tree.resample_nearest(pp1, np.flipud(np.moveaxis(tp.values,0,-1)), pp, radius_of_influence=1000),-1,0)

    #plt.imshow(tp1)
    #plt.colorbar()

    #plt.figure()
    #plt.scatter(sm.flatten())
