import os 

jobIDs = [16885197, 16885196, 16885185, 16885198, 16885199, 16885200, 16885201, 16885202, 16885203, 16885204, 16885205, 16885206, 16885207, 16885208, 16885209, 16885210, 16885211, 16885216, 16885219, 16885220, 16885226, 16885228, 16885230, 16885231, 16885232, 16885233, 16885235, 16885236, 16885237, 16885239, 16882527, 16882531, 16882544, 16882553, 16882558, 16882579, 16882581, 16882583, 16882586, 16882588, 16882597, 16882599, 16882600, 16882601, 16882608]

for i in jobIDs:
	os.system('scancel %d'%i)