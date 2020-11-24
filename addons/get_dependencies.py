#!/usr/bin/env python3
# coding: utf-8

import sys, os
import json
import argparse

import requests
import hashlib
import tarfile
import zipfile

def param_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--pkg', dest="pkg", required=True, help='Available packages', choices=PACKAGES)
    parser.add_argument('--arch', dest="arch", required=True, choices=ARCHS, help='Current arch')
    parser.add_argument('--checksum', dest="checksum", default=True, help='Check file integrity after downloading')
    return parser


PACKAGES = ("coin-clp", )
ARCHS = ("linux", "win64", "osx")

def main():
    parser = param_parser()
    args = parser.parse_args()
    json_packages = os.path.join(os.curdir, 'etc')
    json_packages = os.path.join(json_packages, 'packages.json')
    
    packages_dict = {}
    with open(json_packages) as fh:
        packages_dict = json.load(fh)

        
    pkg_dict = packages_dict[args.pkg][args.arch]

    print("Fetching package:", end="")
    print("- %s (%s)" % (args.pkg, args.arch))
    print("Downaling from: %s" % pkg_dict['url'])
    r = requests.get(pkg_dict['url'], allow_redirects=True)

    print("Package cheksum(sha256)", end=" ")
    hash_strn = hashlib.sha256(r.content).hexdigest()
    assert pkg_dict['sha256'] == hash_strn
    print("Ok!")

    if not os.path.exists(args.pkg):
        print("Creating package folder")
        os.mkdir(args.pkg)
    
    fname = args.pkg + ".compress"
    with open(fname, 'wb') as fh:
        fh.write(r.content)

    if args.arch == "win64":
        archiver = zipfile.ZipFile(fname, 'r')
    else:
        archiver = tarfile.open(fname, "r:gz")
 
    print("Extracting package in %s... " % args.pkg, end=" ")
    archiver.extractall(path=args.pkg)
    archiver.close()
    os.remove(fname)
    print("Ok!")
    print("Dependency retrived correctly :-)\n")

if __name__ == "__main__":
    main()