import os
import numpy as np
import scipy
import brainspace.mesh
import brainspace.plotting
import surfplot
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import pyvirtualdisplay

from . import datasets, transform, stats


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
    cbar_kwargs={},
    nan_color=(0.85, 0.85, 0.85, 1),
    as_outline=False,
    embed_nb=True,
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
        kwargs passed to plotter
    **cbar_kwargs
        kwargs passed to colorbar plotter
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
    else:
        screenshot = False
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
            fig = plot_colorbar(vrange[0], vrange[1], c, **cbar_kwargs)
            if filename:
                fig.savefig('.'.join(filename.split('.')[:-1])+'_cbar.svg', bbox_inches='tight', pad_inches=0, dpi=1200)
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

def reg_plot(
    X, Y, parcellation_name, n_perm=1000,
    xlabel='', ylabel='', metrics=['pearson', 'cosine'],
    ax=None,
):
    """
    Performs spin test between X and Y using `tests` and 
    reports them on a regression plot.

    Parameters
    ----------
    X, Y: (np.ndarray)
        arrays to be compared
    parcellation_name: (str)
    n_perm: (int)
    xlabel, ylabel: (str)
    metrics: (list)
        - 'pearson': Pearson correlation
        - 'cosine': Cosine similarity
    ax: (`matplotlib.axes.Axes`)
        if None, a new figure is created
    
    Returns
    -------
    ax: (`matplotlib.axes.Axes`)
    """
    # stats
    statistics = {}
    ps = {}
    if 'pearson' in metrics:
        statistics['pearson'], ps['pearson'], _ = stats.spin_test_parcellated(
            X, Y, parcellation_name, n_perm=n_perm,
            method='pearson'
        )
        print("Correlation coefficient:", statistics['pearson'].iloc[0, 0], "; p-vlaue:", ps['pearson'].iloc[0, 0])
    if 'cosine' in metrics:
        statistics['cosine'], ps['cosine'], _ = stats.spin_test_parcellated(
            X, Y, parcellation_name, n_perm=n_perm, 
            method=lambda a, b: 1 - scipy.spatial.distance.cosine(a, b)
        )
        print("Cosine similarity:", statistics['cosine'].iloc[0, 0], "; p-vlaue:", ps['cosine'].iloc[0, 0])
    # plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(3.5,2.8))
    sns.regplot(
        x=X, 
        y=Y, 
        scatter_kws=dict(color='#44546A', alpha=1, s=10), 
        line_kws=dict(color='#44546A'), 
        ax=ax
    )
    text = ''
    if 'pearson' in metrics:
        text += f"r = {statistics['pearson'].iloc[0,0]:.3f} "
        if ps['pearson'].iloc[0,0] > 0.001:
            text += f"(p = {ps['pearson'].iloc[0,0]:.3f})"
        else:
            text += '(p < 0.001)'
        if ps['pearson'].iloc[0,0] < 0.05:
            text+=r'$^*$'
        text += '\n'
    if 'cosine' in metrics:
        text += f"cos = {statistics['cosine'].iloc[0,0]:.3f} "
        if ps['cosine'].iloc[0,0] > 0.001:
            text += f"(p = {ps['cosine'].iloc[0,0]:.3f})"
        else:
            text += '(p < 0.001)'
        if ps['cosine'].iloc[0,0] < 0.05:
            text+=r'$^*$'
    text = text.replace('p', r'p$_{spin}$')
    text_x = ax.get_xlim()[0]+(ax.get_xlim()[1]-ax.get_xlim()[0])*0.05
    text_y = ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])*1
    ax.text(text_x, text_y, text,
            color='black',size=16,
            multialignment='left')
    ax.set_xlabel(xlabel, fontdict=dict(fontsize=20))
    ax.set_ylabel(ylabel, fontdict=dict(fontsize=20))
    sns.despine()
    return ax

def plot_icc(icc, parcellation_name):
    """
    Plots ICC distribution and surface map

    Parameters
    ----------
    icc: (pd.Series)
        intraclass correlation coefficients across parcels
        Shape: (n_parcels,)
    parcellation_name: (str)
    """
    fig, ax = plt.subplots(figsize=(4, 0.5))
    sns.swarmplot(
        x=icc,
        s=3, color='#44546A',
        ax=ax
    )
    sns.boxplot(
        x=icc,
        showfliers=False,
        showcaps=False, width=0.1,
        boxprops={"facecolor": (1, 1, 1, 0.75)},
        ax=ax)
    ax.set_yticks([])
    ax.set_xlim([-1, 1])
    ax.set_xlabel('Intraclass correlation coefficient')
    sns.despine()
    plt.setp(ax.collections, zorder=0, label="")
    
    return plot_surface(
        transform.deparcellate_surf(icc, parcellation_name, concat=True, space='fsaverage'), 
        'fsaverage', mesh_kind='semi-inflated',
        cmap='magma', cbar=True, vrange=(icc.min().round(2), icc.max().round(2)),
        cbar_kwargs=dict(figsize=(2,2)),
        layout_style='row'
    )

def plot_icc_by_age(icc_by_age, parcellation_name):
    """
    Plots ICC distribution (in all, versus younger and older age groups)
    as well as surface map of ICC across all subjects

    Parameters
    ----------
    icc_by_age: (pd.DataFrame)
        intraclass correlation coefficients across parcels and age groups
        Shape: (n_parcels, 3)
        Columns: 'all', 'younger', 'older'
    parcellation_name: (str)
    """
    t, p = scipy.stats.ttest_rel(icc_by_age['younger'], icc_by_age['older'])
    print(f'T ={t}, p = {p}')
    plot_data = icc_by_age.unstack().reset_index()
    plot_data['level_0'] = plot_data['level_0'].str.title()
    fig, ax = plt.subplots(figsize=(4, 0.85))
    sns.boxplot(
        data=plot_data,
        x=0,
        hue='level_0',
        showfliers=False,
        showcaps=False, width=0.6,
        palette=['lightgrey', 'darkgrey', 'grey'],
        ax=ax,
        legend=True
    )
    ax.legend(loc='center left', fontsize=10, frameon=False)
    ax.set_yticks([])
    ax.set_xlim([-1, 1])
    ax.set_xlabel('Intraclass correlation coefficient')
    sns.despine()
    plt.setp(ax.collections, zorder=0, label="")
    if p < 0.05:
        ax.text(plot_data[0].max()+0.1, -0.125, '|*', va='center', ha='left', color='grey', fontsize=20)
    return plot_surface(
        transform.deparcellate_surf(icc_by_age['all'], parcellation_name, concat=True, space='fsaverage'), 
        'fsaverage', mesh_kind='semi-inflated',
        cmap='magma', cbar=True, vrange=(icc_by_age['all'].min().round(2), icc_by_age['all'].max().round(2)),
        cbar_kwargs=dict(figsize=(2,2)),
        layout_style='row'
    )