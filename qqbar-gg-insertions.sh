# c_1 insertion between leg 1 and 2
./Contract k_[103,3]*k_[104,4]*t_[201,101,1]*t_[201,2,102]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 1 and 2
./Contract k_[103,3]*k_[104,4]*t_[201,101,1]*t_[201,2,102]*d_[103,104,202]*t_[202,102,101]

# c_3 insertion between leg 1 and 2
./Contract k_[103,3]*k_[104,4]*t_[201,101,1]*t_[201,2,102]*c_[0,1]*f_[103,104,202]*t_[202,102,101]

# c_1 insertion between leg 1 and 3
./Contract c_[0,1]*k_[102,2]*k_[104,4]*f_[3,201,103]*t_[201,101,1]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 1 and 3
./Contract c_[0,1]*k_[102,2]*k_[104,4]*f_[3,201,103]*t_[201,101,1]*d_[103,104,202]*t_[202,102,101]

# c_3 insertion between leg 1 and 3
./Contract c_[0,1]*k_[102,2]*k_[104,4]*f_[3,201,103]*t_[201,101,1]*c_[0,1]*f_[103,104,202]*t_[202,102,101]

# c_1 insertion between leg 1 and 4
./Contract c_[0,1]*k_[102,2]*k_[103,3]*f_[4,201,104]*t_[201,101,1]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 1 and 4
./Contract c_[0,1]*k_[102,2]*k_[103,3]*f_[4,201,104]*t_[201,101,1]*d_[103,104,202]*t_[202,102,101]

# c_3 insertion between leg 1 and 4
./Contract c_[0,1]*k_[102,2]*k_[103,3]*f_[4,201,104]*t_[201,101,1]*c_[0,1]*f_[103,104,202]*t_[202,102,101]

# c_1 insertion between leg 3 and 4
./Contract c_[-1,0]*k_[101,1]*k_[102,2]*f_[3,201,103]*f_[4,201,104]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 3 and 4
./Contract c_[-1,0]*k_[101,1]*k_[102,2]*f_[3,201,103]*f_[4,201,104]*d_[103,104,202]*t_[202,102,101]

# c_3 insertion between leg 3 and 4
./Contract c_[-1,0]*k_[101,1]*k_[102,2]*f_[3,201,103]*f_[4,201,104]*c_[0,1]*f_[103,104,202]*t_[202,102,101]
