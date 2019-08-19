import os
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from glob import glob
from scipy.spatial import Voronoi
from matplotlib.collections import LineCollection
from ..definitions import lipidmolecules
from .. import log

LOGGER = log.LOGGER

def _create_images(points: np.array, box: np.array) -> np.array:
    ''' A set of points is copied eightfold using dimensions in box '''
    outp_ar = []
    if not np.all(box != 0):
        raise ValueError("Box edge length must not be 0")
    img_vectors = [
        (0,  1),
        (1,  0),
        (1,  1),
        (0, -1),
        (-1, 0),
        (-1, -1),
        (1, -1),
        (-1, 1),
        (0, 0),
        ]
    for vec in img_vectors:
        newb = vec * box
        newp = points + newb
        outp_ar.append(newp)
    return np.array([p for ar in outp_ar for p in ar])

def _colorize(regions, ax, fillparameter="shape"):
    ''' Fill Voronoi regions(polygons) depending on fillparameter '''
    if fillparameter == "shape":
        colordict = {
            3:"xkcd:brown",
            4:"xkcd:crimson",
            5:"xkcd:coral",
            6:"xkcd:cyan",
            7:"xkcd:aquamarine",
            8:"xkcd:lime",
            9:"xkcd:green",
            10:"xkcd:yellowgreen",
            11:"xkcd:yellow",
            12:"xkcd:wheat",
            13:"xkcd:brown",
            14:"xkcd:chocolate",
            15:"xkcd:black"
        }
        for reg in regions:
            n_edges = reg.shape[0]
            ax.fill(*zip(*reg), color=colordict[n_edges], alpha=0.8)
    else:
        raise NotImplementedError("Available fillparameters are: shape")
    return ax

def point_in_box(point: np.array, box: np.array) -> bool:
    ''' Return True if point lies in box with vectors starting from 0 '''
    if np.all(point >= 0) and np.all(point <= box):
        return True
    else:
        return False

def plot_periodic_voro(points, box, colorfill="shape", plot_points=False, **kw):
    ''' Create Voronoi plot for a set of inputpoints
        -points must be 2D
        -box must be 2D
        The voronoi tiles can be filled using different schemes, by now only coloring
        regarding number of edges is implemented only
    '''

    if isinstance(points, tuple):
        naxs = len(points)
        fig, axs= plt.subplots(ncols=naxs, figsize=(4*naxs, 4) )
    else:
        naxs = 1
        fig, axs= plt.subplots(ncols=1, figsize=(4, 4) )
        axs = [axs]
        points = (points, )

    line_colors = kw.get('line_colors', 'k')
    line_width = kw.get('line_width', 1.0)
    line_alpha = kw.get('line_alpha', 1.0)

    for colindex, p in enumerate(points):
        periodic_points = _create_images(p, box)

        # Get new set of points that all are within box
        points_in_box = []
        for point in periodic_points:
            if point_in_box(point, box):
                points_in_box.append(point)
        points_in_box = np.array(points_in_box)

        # Add points that are inside box
        if plot_points:
            axs[colindex].plot(points_in_box[:,0], points_in_box[:,1], 'k.')

        vor = Voronoi(periodic_points)

        added_points_region = []
        segments_in_box = []
        polygons = []
        # loop over all ridges and point pairs of voronoi graph
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            pts = vor.points[pointidx]
            if np.all(simplex >= 0): # Only take into account finite ridges (-1 means second point of ridge outside)

                # Only draw ridges if at least one point of ridge pair lies within box
                if np.any(np.all(pts >= 0, axis=1)) and np.any(np.all(pts <= box, axis=1)):
                    segment = vor.vertices[simplex]

                    # Find point that lies within box and save its polygon
                    for i, point in enumerate(pts):
                        if point in added_points_region:
                            continue
                        if point_in_box(point, box):
                            poly = vor.vertices[vor.regions[vor.point_region[pointidx[i]]]]
                            polygons.append(poly)
                            segments_in_box.append(segment)

        axs[colindex].add_collection(LineCollection(segments_in_box,
                                         colors=line_colors,
                                         lw=line_width,
                                         alpha=line_alpha,
                                         linestyle='solid'))
        axs[colindex].grid(b=False)
        axs[colindex].axis("off")
        plt.tight_layout()
        if colorfill:
            ax = _colorize(polygons, axs[colindex], fillparameter=colorfill)
    return fig, axs


def plot_voro_on_structure(mda_atoms, selstr, plot_points=False, output_filename="voronoi.png"):
    ''' '''
    colordict = {
        "DPPC":"Red",
        "DUPC":"Black",
        "CHL1":"Yellow",
        "ERG":"Green",
    }
    points = []
    if not isinstance(selstr, list):
        selstr = [selstr]
    for sel in selstr:
        points.append(mda_atoms.select_atoms(sel).positions[:, 0:2])
    box = mda_atoms.dimensions[:2]
    fig, axs = plot_periodic_voro(tuple(points), box, plot_points=False, colorfill="shape")
    if plot_points:
        resnames = set(mda_atoms.resnames) - set(lipidmolecules.SOLVENTS)
        for i, sel in enumerate(selstr):
            for resn in resnames:
                points = mda_atoms.select_atoms(sel+" and resname {}".format(resn)).positions[:, 0:2]
                periodic_points = _create_images(points, box)

                # Get new set of points that all are within box
                points_in_box = []
                for point in periodic_points:
                    if point_in_box(point, box):
                        points_in_box.append(point)
                points_in_box = np.array(points_in_box)
                axs[i].plot(points_in_box[:,0], points_in_box[:,1], '.', color=colordict[resn], markeredgecolor="black")
    fig.savefig(output_filename)
    plt.close()

def make_voro_movie(systeminfo, video_name="voro_movie.mpg", dt=None, lvl="head", leaflet=0, delete_imgs=True, plot_points=True,
    imgdir="movie", overwrite=False,
    ):
    ''' '''
    lvls = {
        "head":"name P O3",
         "tail2":"(resname DPPC and name C12 C22) or (resname CHL1 and name O3)",
         "tail4":"(resname DPPC and name C14 C24) or (resname CHL1 and name O3)",
         "tail6":"(resname DPPC and name C16 C26) or (resname CHL1 and name O3)",
         "tail8":"(resname DPPC and name C18 C28) or (resname CHL1 and name O3)",
        "tail10":"(resname DPPC and name C110 C210) or (resname CHL1 and name O3)",
        "tail12":"(resname DPPC and name C112 C212) or (resname CHL1 and name O3)",
        "tail14":"(resname DPPC and name C114 C214) or (resname CHL1 and name O3)",
    }
    if lvl not in lvls.keys():
        atomselect_str = lvl
    else:
        atomselect_str = lvls[lvl]

    IMG_TMPDIR = imgdir

    if isinstance(leaflet, int):
        leaflet = (leaflet,)
    elif isinstance(leaflet, list):
        leaflet = tuple(leaflet)

    selstr = []
    for leaf in leaflet:
        resids =  []
        if dt is None or dt < systeminfo.dt:
            dt = systeminfo.dt
        for i in systeminfo.MOLRANGE:
            if systeminfo.res_to_leaflet[i] == leaf: resids.append(i)

        ## 2 Create selection
        selstr.append('{} and resid {}'.format(atomselect_str, ' '.join([str(i) for i in resids])))

    ## 1 Create Folder
    os.makedirs(IMG_TMPDIR, exist_ok=True)
    #try:
    #    os.makedirs(IMG_TMPDIR, exist_ok=overwrite)
    #except FileExistsError:
    #    cnt = 0
    #    while os.path.exists(IMG_TMPDIR+str(cnt)):
    #        cnt += 1
    #    IMG_TMPDIR = IMG_TMPDIR+str(cnt)
    #    os.makedirs(IMG_TMPDIR, exist_ok=False)

    ## Print logfile: which information are stored in video??
    with open(IMG_TMPDIR+"/info.log", "w") as f:
        print(systeminfo.info(), file=f)
        print("Video settings:", file=f)
        print("dt:", dt, file=f)
        print("lvl:", atomselect_str, file=f)
        print("leaf:", leaflet, file=f)

    ## 3 Loop over trajectory with dt and create voro in tmp folder
    len_traj = len(systeminfo.universe.trajectory)
    Nzeros = str(int(np.ceil(np.log10(len_traj))))
    picture_filename_template =  "./{}/{:0Nd}_voro.png".replace("N", Nzeros)
    for t in range(len_traj):
        picture_filename = picture_filename_template.format(IMG_TMPDIR, t)
        time = systeminfo.universe.trajectory[t].time
        LOGGER.info("At time %s", time)
        if time % dt != 0 or systeminfo.t_start > time:
            continue
        elif systeminfo.t_end < time:
            break
        if os.path.isfile(picture_filename) and not overwrite:
            continue
        plot_voro_on_structure(systeminfo.universe.atoms,
            selstr,
            plot_points=plot_points,
            output_filename=picture_filename,
            )

    ## 4 Gather voropicturefiles and create movie
    cmd = ["ffmpeg", "-framerate", "2",
        "-i", "{}/%0{}d_voro.png".format(IMG_TMPDIR, Nzeros),
        "-c:v", "mpeg2video",
        "-pix_fmt", "yuv420p",
        "{}/{}".format(IMG_TMPDIR, video_name)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    proc.wait()
    proc.stdout.close()
    proc.stderr.close()
    print(out.decode())
    print(err.decode())

    ## OPTIONALLY: Delete img files
    if delete_imgs:
        for png in glob(IMG_TMPDIR+"/*.png"):
            os.remove(png)
