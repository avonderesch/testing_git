
n_elements = 10   ; Number of elements in the array
min_temp = -9.0  ; Minimum temperature value
max_temp = 15.0   ; Maximum temperature value

min_prec = 0  ; Minimum temperature value
max_prec = 15.0   ; Maximum temperature value

; Marit and Matthias
<<<<<< HEAD
; hello from alex and hello too from Lander and also from the rest in
; the room and whole Sion
=======
; hello from alex and hello too from Lander !!
>>>>>>> 66c375c91bf85f4e351a1dbf7fad9b6f935d619e
; Generate random numbers between 0 and 1
random_numbers = RANDOMU(seed, n_elements)

; Scale and shift random numbers to represent temperature values
temperature_values = min_temp + random_numbers * (max_temp - min_temp)
precipitation_values = min_prec + random_numbers * (max_prec - min_prec)

snow=0
; hello from Lander
; hello from Friedrich
; And too
for i=0, 9 do begin

print, i

temp=temperature_values(i)
prec=precipitation_values(i)

if temp lt 0 then snow=prec else snow=0
; hello from Lander
snow=snow+snow

endfor

print, temperature_values
print, precipitation_values

print, "this is the amount of snow:"
print, snow

end 


