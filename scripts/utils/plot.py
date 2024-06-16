import os
import numpy as np
import brainspace.mesh
import brainspace.plotting
import surfplot
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import pyvirtualdisplay

from . import datasets


def plot_surface(
    surface_data,
    mesh,
    itype=None,
    mesh_kind="inflated",
    filename=None,
    layout_style="row",
    cmap="viridis",
    vrange=None,
    cbar=False,
    nan_color=(0.85, 0.85, 0.85, 1),
    as_outline=False,
    plotter="brainspace",
    **plotter_kwargs,
):
    """
    Plots `surface_data` on `mesh` using brainspace

    Parameters
    ----------
    surface_data: (np.ndarray)
        for brainspace must be 1D, for surfplot can be nD (m, vert) (will be plotted as layers)
    mesh: (str | dict)
        - fsaverage
        - fsaverage5
        - fsLR (32k)
        - dict of path to meshes for 'L' and 'R'
    itype: (str | None)
        mesh file type. For .gii enter None. For freesurfer files enter 'fs'
    mesh_kind: (str)
        - 'orig'
        - 'inflated'
        - 'sphere'
    filename: (Pathlike str)
    cmap: (str | list of str)
        cmap must be a list if plotter is surfplot and multiple arrays are provided
    layout_style: (str)
        - row
        - grid
    vrange: (tuple | None)
    cbar: (bool)
    nan_color: (tuple)
    as_outline: (bool or list of bool)
        only used with surfplot
    plotter: {'brainspace', 'surfplot'}
    **plotter_kwargs
    """
    if isinstance(cmap, str):
        cmap = [cmap]
    if isinstance(as_outline, bool):
        as_outline = [as_outline]
    # load surface mesh
    if isinstance(mesh, str):
        downsampled = mesh == "fsaverage5"
        if mesh.startswith("fsaverage"):
            space = "fsaverage"
        else:
            space = mesh
        mesh = datasets.load_mesh_paths(
            kind=mesh_kind, space=space, downsampled=downsampled
        )
    for fs_suffix in [
        ".pial",
        ".midthickness",
        ".white",
        ".inflated",
        ".pial_semi_inflated",
    ]:
        if mesh["L"].endswith(fs_suffix):
            itype = "fs"
    if not os.path.exists(mesh["L"]):
        raise ValueError("Mesh not found")
    lh_surf = brainspace.mesh.mesh_io.read_surface(mesh["L"], itype=itype)
    rh_surf = brainspace.mesh.mesh_io.read_surface(mesh["R"], itype=itype)
    # configurations
    if filename:
        screenshot = True
        embed_nb = False
        filename += ".png"
    else:
        screenshot = False
        embed_nb = True
    if layout_style == "row":
        size = (1600, 400)
        zoom = 1.2
    else:
        size = (900, 500)
        zoom = 1.8
    if vrange is None:
        vrange = (np.nanmin(surface_data), np.nanmax(surface_data))
    elif vrange == "sym":
        vmin = min(np.nanmin(surface_data), -np.nanmax(surface_data))
        vrange = (vmin, -vmin)
    if cbar:
        for c in cmap:
            plot_colorbar(vrange[0], vrange[1], c)
    # create virtual display for plotting in remote servers
    disp = pyvirtualdisplay.Display(visible=False)
    disp.start()
    if plotter == "brainspace":
        return brainspace.plotting.surface_plotting.plot_hemispheres(
            lh_surf,
            rh_surf,
            surface_data,
            layout_style=layout_style,
            cmap=cmap[0],
            color_range=vrange,
            # TODO: change size and zoom based on layout
            size=size,
            zoom=zoom,
            interactive=False,
            embed_nb=embed_nb,
            screenshot=screenshot,
            filename=filename,
            transparent_bg=True,
            nan_color=nan_color,
            **plotter_kwargs,
        )
    elif plotter == "surfplot":
        # initialize the plot
        p = surfplot.Plot(
            surf_lh=lh_surf,
            surf_rh=rh_surf,
            layout=layout_style,
            size=size,
            zoom=zoom,
            brightness=1,
            mirror_views=True,
        )
        # add the layers
        for i in range(surface_data.shape[0]):
            p.add_layer(
                surface_data[i, :], cmap=cmap[i], as_outline=as_outline[i], cbar=False
            )
        # build
        p.build()
        return p


def plot_colorbar(
    vmin, vmax, cmap=None, bins=None, orientation="vertical", figsize=None
):
    """
    Plots a colorbar
    Parameters
    ---------
    vmin, vmax: (float)
    cmap: (str or `matplotlib.colors.Colormap`)
    bins: (int)
        if specified will plot a categorical cmap
    orientation: (str)
        - 'vertical'
        - 'horizontal'
    figsize: (tuple)

    Returns
    -------
    fig: (`matplotlib.figure.Figure`)
    """
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        np.linspace(vmin, vmax, 100).reshape(10, 10), cmap=plt.cm.get_cmap(cmap, bins)
    )
    fig.gca().set_visible(False)
    divider = make_axes_locatable(ax)
    if orientation == "horizontal":
        cax = divider.append_axes("bottom", size="10%", pad=0.05)
    else:
        cax = divider.append_axes("left", size="10%", pad=0.05)
    fig.colorbar(im, cax=cax, ticks=np.array([vmin, vmax]), orientation=orientation)
    cax.yaxis.tick_left()
    cax.xaxis.tick_bottom()
    return fig


def plot_dynamic_fc(window_fcs, output_path):
    """
    Creates a movie of dynamic FC

    Parameters
    ----------
    windwo_fcs: (np.ndarray)
        dynamic FC with shape (nodes, nodes, windows)
    output_path: (str)
        path to a .gif or .mp4 file
    """
    import imageio

    tmp_img_path = f"/tmp/tmp_{output_path}.png"
    images = []
    for window in range(window_fcs.shape[2]):
        fig, ax = plt.subplots(1)
        sns.heatmap(window_fcs[:, :, window], vmin=-1, vmax=1, ax=ax)
        ax.set_title(f"Window: {window+1}")
        fig.savefig(tmp_img_path)
        plt.close(fig)
        images.append(imageio.imread(tmp_img_path))
    imageio.mimsave(output_path, images, fps=5)
