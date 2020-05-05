#!/usr/bin/env python3
import warnings
import gzip
import bz2
from zipfile import ZipFile
import io


def parse_vars(file_handle):
    skip_lines = 0
    for line in file_handle:
        line = line.rstrip()
        if line.startswith("#"):
            skip_lines = skip_lines + 1
        else:
            break
    return skip_lines


def determine_header_structure(file):
    f = open(file)
    line = f.readline()
    header_exists = None
    if line.startswith("[Header]") or line.startswith("[Data]"):
        header_exists = True
    else:
        header_exists = False
    f.close()
    return header_exists


def parse_header(file):
    does_header_exist = determine_header_structure(file)
    if does_header_exist is True:
        with open(file, 'r') as input_file:
            header_dict = {}
            count = 1
            header_skip = 0
            for line in input_file:
                if line.startswith("[Data]"):
                    header_skip = count
                    break
                else:
                    if not line.startswith("[Header]"):
                        line = line.rstrip("\n")
                        line_array = line.split("\t")
                        if line_array[0] == 'Content':
                            header_dict.update({line_array[0]: line_array[2]})
                        else:
                            header_dict.update({line_array[0]: line_array[1]})
                    count = count + 1
        input_file.close()
    else:
        print("[Header]...[Data]... file structure not detected. Assuming first line contains data column names")
        header_skip = 0
        header_dict = {}
    if does_header_exist is True and header_skip == 0:
        print("[Header]...[Data]... file structure not detected. Assuming first line contains data column names")
    else:
        pass
    return header_skip, header_dict


# Information for determining the type of file compression
magic_dict = {
    "\x1f\x08\x08": "gzip",
    "\x1f\x8b\x08": "gzip",
    "\x42\x5a\x68": "bz2",
    "\x50\x4b\x03\x04": "zip",
    "\x53\x70\x5f": "tar",
    "\x75\x73\x74\x61\x62": "tar"
    }

max_len = max(len(x) for x in magic_dict)


def get_compression_type(path):
    with open(path, "r", errors="ignore", encoding="utf-8") as file_path:
        file_start = file_path.read()
        for magic, filetype in magic_dict.items():
            if file_start.startswith(magic):
                return filetype
        return "no match"


def get_uncompressed_file(path):
    # if an uncompressed text file
    if path.endswith((".txt", ".csv")):
        filetype = "native_uncompressed"
    # if a not a text or csv file (potentially compressed)
    else:
        compression_type = get_compression_type(path)
        # compressed and has an expected file extension
        if path.endswith((".gzip", ".gz", ".bz2", ".zip", ".xz")):
            # compressed but with tar
            if path.endswith(".tar.gz"):
                warnings.warn("tar compressed files are not allowed - un-tar and compress using gzip, zip, bzip2,or xz",
                              stacklevel=4)
                filetype = "skip_file"
            # file has an acceptable file extension, test for compression type
            else:
                compression_type = get_compression_type(path)
                if compression_type == "gzip" and path.endswith((".gz", ".gzip")):
                    filetype = "gzip"
                elif compression_type == "bz2" and path.endswith((".bz2", ".bzip2")):
                    filetype = "bz2"
                elif compression_type == "zip" and path.endswith(".zip"):
                    filetype = "zip"
                elif compression_type == "tar":
                    warnings.warn("File " + path + " may be compressed with tar - skipping")
                    filetype = "skip_file"
                elif compression_type == "xz" and path.endswith(".xz"):
                    filetype = "xz"
                else:
                    filetype = compression_type
                    warnings.warn("File extension may not match determined compression type " + compression_type,
                                  stacklevel=4)
        # potentially compressed but does not have an expected extension
        else:
            warnings.warn(path + " does not have a recognized file extension", stacklevel=4)
            # compression type might be identified by magic number
            if compression_type != "no match":
                warnings.warn("File " + path + " may be in " + compression_type + " format - uncompressing",
                              stacklevel=4)
                filetype = compression_type
            # compression type cannot be identified by magic number; may not be compressed, may be corrupted, etc.
            else:
                warnings.warn("File type for " + path + " cannot be determined - skipping", stacklevel=4)
                filetype = "skip_file"
    return filetype


def uncompressing(path):
    filetype = get_uncompressed_file(path)
    if filetype == "no match" or filetype == "native_uncompressed":
        # we want to actually open the file here
        with open(path, "r", errors="ignore", encoding="utf-8") as uncompressed_handle:
            return parse_vars(uncompressed_handle)
    elif filetype == "gzip":
        with gzip.open(path, "rt", errors="ignore") as gzipped_handle:
            return parse_vars(gzipped_handle)
    elif filetype == 'bz2':
        with bz2.open(path, "rt", errors="ignore") as bzipped_handle:
            return parse_vars(bzipped_handle)
    elif filetype == "zip":
        with ZipFile(path, mode='r', allowZip64=True) as zipped_file:
            zip_list = ZipFile.namelist(zipped_file)
            if len(zip_list) != 1:
                warnings.warn("Zipped archive " + path + " contains more than one file.")
            else:
                with zipped_file.open(zip_list[0], mode='r') as zipped_handle:
                    items_file = io.TextIOWrapper(zipped_handle, encoding='UTF-8', newline='')
                    return parse_vars(items_file)
    else:
        pass


def check_header_data_congruency(header_dict, calc_num_samples, calc_num_snps):
    warning_message_list = []
    if int(header_dict['Num Samples']) != calc_num_samples:
        message = "Inconsistent number of samples in input file based on 'Num Samples' in header"
        warning_message_list.append(message)
    else:
        pass
    if int(header_dict['Num SNPs']) != calc_num_snps:
        message = "Inconsistent number of SNPs in input file based on 'Num SNPs' in header"
        warning_message_list.append(message)
    else:
        pass
    # print warnings
    for warning_entry in warning_message_list:
        print(warning_entry)
    return warning_message_list
