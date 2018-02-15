kd_sort = function(x, ...) UseMethod("kd_sort")

kd_sort.matrix = function(x)
{
  y = matrix_to_tuples(x)
  kd_sort_(y, inplace = TRUE)
  return(tuples_to_matrix(y))
}

kd_sort.arrayvec = function(x, inplace = FALSE)
{
  return(kd_sort_(x, inplace = inplace))
}

kd_lower_bound = function(x, v) UseMethod("kd_lower_bound")

kd_lower_bound.matrix = function(x, v)
{
  y = matrix_to_tuples(x)
  return(kd_lower_bound_(y, v))
}

kd_lower_bound.arrayvec = function(x, v)
{
  return(kd_lower_bound_(x, v))
}

kd_upper_bound = function(x, v) UseMethod("kd_upper_bound")

kd_upper_bound.matrix = function(x, v)
{
  y = matrix_to_tuples(x)
  return(kd_upper_bound_(y, v))
}

kd_upper_bound.arrayvec = function(x, c)
{
  return(kd_upper_bound_(x, v))
}

kd_range_query = function(x, l, u) UseMethod("kd_range_query")

kd_range_query.matrix = function(x, l, u)
{
  y = matrix_to_tuples(x)
  z = kd_range_query_(y, l, u)
  return(tuples_to_matrix(z))
}

kd_range_query.arrayvec = function(x, l, u)
{
  return(kd_range_query_(x, l, u))
}
