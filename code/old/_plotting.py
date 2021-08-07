def get_overlay(roi_ind, value, hemi):

  is_rh = roi_ind > 200

  if hemi == "right":

    inds = np.array(roi_ind[is_rh] - 200)
    vals = value[is_rh]
    overlay = np.zeros(len(schaefer_rh[0]))
    template = schaefer_rh[0]

  elif hemi == "left":

    inds = np.array(roi_ind[~is_rh])
    vals = value[~is_rh]
    overlay = np.zeros(len(schaefer_lh[0]))
    template = schaefer_lh[0]

  for i in range(len(inds)):

    roi_i = inds[i]
    val_i = vals.iloc[np.nonzero(inds == roi_i)]
    overlay[template == roi_i] = val_i

  return overlay



def plot_surf_roi_montage(
  roi_map, title, fs_surf_mesh = 'pial', fs_bg_map = 'sulc', bg_on_data = False, darkness = 1, cmap = 'magma',
  ):

  hemis = ("left", "right")
  views = ("lateral", "medial")

  fig, axs = plt.subplots(1, 4, subplot_kw = dict(projection = "3d"), figsize = plt.figaspect(1/4))

  for hemi_i in range(len(hemis)):
    for view_i in range(len(views)):
      # hemi_i = 1
      # view_i = 1

      plotting.plot_surf_roi(

        fsaverage[fs_surf_mesh + '_' + hemis[hemi_i]],
        roi_map = roi_map[:, hemi_i],
        hemi = hemis[hemi_i],
        view = views[view_i],
        bg_map = fsaverage[fs_bg_map + '_' + hemis[hemi_i]],
        bg_on_data = bg_on_data,
        darkness = darkness,
        cmap = cmap,
        colorbar = view_i == 1 & hemi_i == 1,
        axes = axs[hemi_i*2 + view_i]

        )

      plt.tight_layout(pad = 0, w_pad = 0, h_pad = 0)

  fig.suptitle(title, fontsize = "xx-large")
  return fig
