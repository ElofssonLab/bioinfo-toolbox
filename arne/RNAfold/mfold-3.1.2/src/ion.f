c     Isolate [Na+] and [Mg++] correction into a separate file.

      real function ion(na_type,olipoly,type,t,na_conc,mg_conc)
      real alog,sqrt,t,na_conc,mg_conc
      character na_type*3,olipoly*7,type*6

      if (na_type.eq.'RNA'.or.na_type.eq.'rna') then
         ion = 0.0
      elseif (olipoly.eq.'polymer'.or.olipoly.eq.'POLYMER') then
         if (type.eq.'affine') then
            ion = -(273.15 + t)/310.15*(0.2 + 0.175*alog(na_conc))
         else
            ion = -(273.15 + t)/310.15*(0.175*alog(na_conc))
         endif
      else
         ion = -(273.15 + t)/310.15*0.114*alog(na_conc + 3.3*sqrt(mg_conc))
      endif
      return
      end

