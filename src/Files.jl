function does_data_file_exist(path_to_data_file::String)::Bool

    if (ispath(path_to_data_file) == true)
        return true
    end

    return false
end