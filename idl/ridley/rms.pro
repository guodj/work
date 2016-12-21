
function rms, array

  result = sqrt(mean(array^2))
  return, result

end
