
pro strsub, name, subout, subin

  name = strmid(name,0,strpos(name,subout))+subin

end
