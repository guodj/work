pro ppos, pos

  xs = float(!d.x_size)
  ys = float(!d.y_size)
  pos = float(!p.clip)
  pos = pos(0:3)
  pos([0,2]) = pos([0,2])/xs
  pos([1,3]) = pos([1,3])/ys

  return

end
