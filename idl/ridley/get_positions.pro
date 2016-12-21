
;------------------------------------------------------------------------------
;
;  Plot the variables subroutines
;
;------------------------------------------------------------------------------


PRO Get_Positions, plottype, placement, numofgra, numofvar, titlenam, screen

    type = plottype

;
; Plot type effected :
;   type_1
;   type_3
;

    if type eq 3 then type = 1

;
; Plot type effected :
;   type_8
;   type_9
;   type_13
;

    if (type eq 9) then type = 8
    if (type eq 13) then type = 8

    case (type) of

      1 : begin

        numdone = 0

        commenu = [' Choose Order of variables : ', 			      $
		   '          '+titlenam,       			      $
                   '   Order rest as displayed   ']

	select = 0

	order = indgen(numofvar+1)*0 -1

	while (select lt numofvar+2) and (numdone le numofvar) do begin

          if n_elements(screen) eq 0 then 			$
	    select = wmenu(commenu, title=0, init=1)		$
	  else begin
	    select = 0
	    while (select le 0) or (select gt n_elements(commenu)-1) do begin
	      print, commenu(0)
	      for j=1,n_elements(commenu)-1 do			$
		print, tostr(j),'. ',commenu(j)
	      print, 'Enter choice:'
	      read, select
	    endwhile
	  endelse

          if select ne numofvar+2 then begin

            if order(select-1) ne -1 then begin

	      commenu(select) = '          '+titlenam(select-1)
	      reorder = where(order gt order(select-1), count)

	      if count gt 0 then begin

		for i=0,count-1 do begin

	          order(reorder(i)) = order(reorder(i)) - 1
                  if order(reorder(i)) ge 9 then 			      $
		    commenu(reorder(i)+1) =                                   $
		      strcompress(string(order(reorder(i))+1),/remove_all)+   $
		      '        '+titlenam(reorder(i)) 			      $
	          else                                                        $
		    commenu(reorder(i)+1) = '0'+                              $
		      strcompress(string(order(reorder(i))+1),/remove_all)+   $
		      '        '+titlenam(reorder(i))

                endfor

	      endif

	      order(select-1) = -1

	      numdone = numdone - 1

	    endif else begin

	      order(select-1) = numdone

              if numdone ge 9 then 					      $
		commenu(select) =                                             $
		  strcompress(string(numdone+1),/remove_all)+ 		      $
		  '        '+titlenam(select-1) 			      $
	      else                                                            $
		commenu(select) = '0'+                                        $
		  strcompress(string(numdone+1),/remove_all)+  		      $
		  '        '+titlenam(select-1)

	      numdone = numdone + 1

            endelse

	  endif

	endwhile

        if numdone le numofvar then begin

          addorder = where(order eq -1, count)

	  if count gt 0 then for i=numdone,numofvar do	 		      $
	    order(addorder(i-numdone)) = i

	endif

        if plottype eq 1 then begin

          for curvar = 0, numofvar do 					      $
	    placement(curvar,0) = where(order eq curvar)
          numofgra = numofvar

        endif else begin

          for curvar = 0, numofvar do                   		      $
	    placement(0,curvar) = where(order eq curvar)
          numofgra = 0

        endelse

      end

      2 : begin

        commenu = [' Choose Variables for graph 01 (4 Max) : ',		      $
		   '                '+titlenam,       			      $
                   '                Next Graph               ']

        numselected = 0
        curgraph    = 0
        i           = 0

        while numselected ne numofvar+1 do begin

          if n_elements(screen) eq 0 then begin
	    placement(curgraph, i) = wmenu(commenu, title=0, init=1)-1
	  endif else begin
	    select = 0
	    while (select le 0) or (select gt n_elements(commenu)-1) do begin
	      print, commenu(0)
	      for j=1,n_elements(commenu)-1 do			$
		print, tostr(j),'. ',commenu(j)
	      print, 'Enter choice:'
	      read, select
	    endwhile
	    placement(curgraph,i) = select-1
	  endelse

	  if (placement(curgraph, i) eq numofvar+1) or (i eq 4) then begin

            if i ne 0 then begin

 	      placement(curgraph, i) = placement(curgraph, i-1)
 	      i                      = 0
	      curgraph               = curgraph + 1
	      commenu(0) = 'Choose Variables for graph '+		      $
		   strcompress(string(curgraph+1),/remove_all)+' (4 Max) :'

	    endif

	  endif else begin

	    commenu(placement(curgraph, i)+1) = 'Graph '+		      $
		  strcompress(string(curgraph+1),/remove_all)+'        '+     $
		  titlenam(placement(curgraph, i))
	    i = i + 1
	    numselected = numselected + 1

	  endelse

        endwhile

        placement(curgraph, i)   = placement(curgraph, i-1)
        placement(curgraph, i+1) = 0
 
        numofgra = curgraph

      end

      4 : begin

        commenu = ['Choose 2 Variables for graph 01 :','             '+titlenam]

        numselected = 0
        curgraph    = 0
        i           = 0

        while numselected ne numofvar+1 do begin

          if n_elements(screen) eq 0 then begin
	    placement(curgraph, i) = wmenu(commenu, title=0, init=1)-1
	  endif else begin
	    select = 0
	    while (select le 0) or (select gt n_elements(commenu)-1) do begin
	      print, commenu(0)
	      for j=1,n_elements(commenu)-1 do			$
		print, tostr(j),'. ',commenu(j)
	      print, 'Enter choice:'
	      read, select
	    endwhile
	    placement(curgraph,i) = select-1
	  endelse

	  if i eq 1 then begin

	    commenu(placement(curgraph, i)+1) = 'Graph '+		      $
		strcompress(string(curgraph+1),/remove_all)+'     '+	      $
		titlenam(placement(curgraph, i))
 	    placement(curgraph, i+1)          = placement(curgraph, i)
            placement(curgraph, i+2)          = 0
	    numselected                       = numselected + 1
 	    i                                 = 0
	    curgraph                          = curgraph + 1
	    commenu(0)                        = 			      $
		'Choose 2 Variables for graph '+			      $
		strcompress(string(curgraph+1),/remove_all)+' :'

	  endif else begin

	    commenu(placement(curgraph, i)+1) = 'Graph '+		      $
		strcompress(string(curgraph+1),/remove_all)+'     '+	      $
		titlenam(placement(curgraph, i))
	    i = i + 1
	    numselected = numselected + 1

	  endelse

        endwhile

        numofgra = curgraph - 1

      end

      8 : begin

        if plottype eq 8 then						      $
          vars = ['Range    ','Azimuth  ','Velocity ','Elevation','Event    ']
        if plottype eq 9 then						      $
          vars = ['Range    ','Azimuth  ','Velocity ','Elevation']
        if plottype eq 13 then						      $
          vars = ['Range    ','Azimuth  ','LOS Vel. ',			      $
		  'Theta Vel','Elevation','Event    ']

        if numofvar eq n_elements(vars)-1 then begin

          commenu = [' Choose variable to represent ' + vars(0) + ' : ',      $
		     '            '+titlenam]

	  num = indgen(numofvar+1)
	  order = indgen(numofvar+1)*0 - 1
          numdone = 0

	  while (numdone le numofvar) do begin

            if n_elements(screen) eq 0 then 			$
	      select = wmenu(commenu, title=0, init=1)		$
	    else begin
	      select = 0
	      while (select le 0) or (select gt n_elements(commenu)-1) do begin
	        print, commenu(0)
	        for j=1,n_elements(commenu)-1 do			$
		  print, tostr(j),'. ',commenu(j)
	        print, 'Enter choice:'
	        read, select
	      endwhile
	    endelse

            next = where(num ne -1, count)
	    commenu(select) = vars(num(next(0)))+' - ' + titlenam(select-1)

            if order(select-1) ne -1 then begin

              num(order(select-1)) = order(select-1)

	    endif else begin

	      numdone = numdone + 1

	    endelse

	    num(next(0)) = -1
	    order(select-1) = next(0)

            next = where(num ne -1, count)

	    if count gt 0 then 						      $
	      commenu(0) = ' Choose variable to represent ' + 		      $
			   vars(num(next(0))) + ' : '

	  endwhile

          for curvar = 0, numofvar do                   		      $
	    placement(0,curvar) = where(order eq curvar)
          numofgra = 0

        endif 

      end

      22 : begin

        commenu = [' Choose Variable to represent X axis : ',		      $
		   '              '+titlenam]

          if n_elements(screen) eq 0 then 			$
	    i = wmenu(commenu, title=0, init=1)	-1		$
	  else begin
	    i = 0
	    while (i le 0) or (i gt n_elements(commenu)-1) do begin
	      print, commenu(0)
	      for j=1,n_elements(commenu)-1 do			$
		print, tostr(j),'. ',commenu(j)
	      print, 'Enter choice:'
	      read, i
	      i = i - 1 
	    endwhile
	  endelse

	numofgra = numofvar-1
	curgraph = 0

	for j=0,numofgra+1 do begin

	  if j ne i then begin

	    placement(curgraph, 0) = i
	    placement(curgraph, 1) = j
	    placement(curgraph, 2) = j
	    placement(curgraph, 3) = 0
	    curgraph = curgraph + 1

	  endif

	endfor

      end

      else : begin

          for curvar = 0, numofvar do                   		      $
	    placement(0,curvar) = curvar
          numofgra = 0
	
      end

;      5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20: begin

;        print, 'in sub'

;      end

      endcase

  RETURN

END
