# c_1 insertion between leg 1 and 2
./Spectrum c_[-1,0]*k_[103,3]*k_[104,4]*f_[1,201,101]*f_[2,201,102]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 1 and 2
./Spectrum c_[-1,0]*k_[103,3]*k_[104,4]*f_[1,201,101]*f_[2,201,102]*d_[101,102,202]*t_[202,103,104]

# c_3 insertion between leg 1 and 2
./Spectrum c_[-1,0]*k_[103,3]*k_[104,4]*f_[1,201,101]*f_[2,201,102]*c_[0,1]*f_[101,102,202]*t_[202,103,104]

# c_1 insertion between leg 1 and 3
./Spectrum c_[0,1]*k_[102,2]*k_[104,4]*f_[1,201,101]*t_[201,3,103]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 1 and 3
./Spectrum c_[0,1]*k_[102,2]*k_[104,4]*f_[1,201,101]*t_[201,3,103]*d_[101,102,202]*t_[202,103,104]

# c_3 insertion between leg 1 and 3
./Spectrum c_[0,1]*k_[102,2]*k_[104,4]*f_[1,201,101]*t_[201,3,103]*c_[0,1]*f_[101,102,202]*t_[202,103,104]

# c_1 insertion between leg 1 and 4
./Spectrum c_[0,1]*k_[102,2]*k_[103,3]*f_[1,201,101]*t_[201,104,4]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 1 and 4
./Spectrum c_[0,1]*k_[102,2]*k_[103,3]*f_[1,201,101]*t_[201,104,4]*d_[101,102,202]*t_[202,103,104]

# c_3 insertion between leg 1 and 4
./Spectrum c_[0,1]*k_[102,2]*k_[103,3]*f_[1,201,101]*t_[201,104,4]*c_[0,1]*f_[101,102,202]*t_[202,103,104]

# c_1 insertion between leg 3 and 4
./Spectrum k_[101,1]*k_[102,2]*t_[201,3,103]*t_[201,104,4]*k_[101,102]*k_[103,104]

# c_2 insertion between leg 3 and 4
./Spectrum k_[101,1]*k_[102,2]*t_[201,3,103]*t_[201,104,4]*d_[101,102,202]*t_[202,103,104]

# c_3 insertion between leg 3 and 4
./Spectrum k_[101,1]*k_[102,2]*t_[201,3,103]*t_[201,104,4]*c_[0,1]*f_[101,102,202]*t_[202,103,104]
