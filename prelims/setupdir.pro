pro setupdir

;does directory structure setup

file_mkdir,'data'
file_mkdir,'dataf'
file_mkdir,'reduc'
file_mkdir,'reduc/expand'
file_mkdir,'reduc/reg'
file_mkdir,'reduc/rsub'
file_mkdir,'reduc/proc'

file_mkdir,'data/raw'
file_mkdir,'data/prep'
file_mkdir,'data/basic'
file_mkdir,'data/cal/dark'
file_mkdir,'data/cal/flat'


end
