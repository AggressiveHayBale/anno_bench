#!/usr/bin/env bash 

fasta=$1
genus=$2 

exec 3<> meta.yaml
    echo "organism:">&3
    echo "  genus_species: '${genus}'">&3
    echo "authors:">&3
    echo "    - author:">&3
    echo "        first_name: 'Firstname'">&3
    echo "        last_name: 'Lastname'">&3
    echo "contact_info:">&3
    echo "        first_name: 'Firstname'">&3
    echo "        last_name: 'Lastname'">&3
    echo "        email: 'Email@address.net'">&3
    echo "        organization: 'Organization'">&3
    echo "        department: 'Department'">&3
    echo "        phone: '301-555-0245'">&3
    echo "        street: '100 Street St'">&3
    echo "        city: 'City'">&3
    echo "        postal_code: '12345'">&3
    echo "        country: 'Country'">&3 
exec 3>&-

exec 3<> input.yaml
    echo "fasta:">&3
    echo "  class: File">&3
    echo "  location: ${fasta}">&3
    echo "submol:">&3
    echo "  class: File">&3
    echo "  location: meta.yaml">&3
    echo "supplemental_data: { class: Directory, location: /pgap/input }">&3
    echo "report_usage: false">&3
    echo "no_internet: true">&3
exec 3>&-