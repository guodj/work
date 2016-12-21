pro get_plot_pos, pos

  pos = fltarr(4)

  pos([0,2]) = float(!p.clip([0,2])) / float(!d.x_size)
  pos([1,3]) = float(!p.clip([1,3])) / float(!d.y_size)

end



