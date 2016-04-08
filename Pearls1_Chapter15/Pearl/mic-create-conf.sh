#!/bin/bash
### Script to create /etc/mpss/mic*.conf
### Currently setup for MPSS 3.2

if [ -d /sys/class/mic ]
then
  for MICN in $( cd /sys/class/mic ; echo mic* )
  do
    MICIPADDR=$(grep ${HOSTNAME}-${MICN} /etc/hosts | awk '{print $1}')
    HOSTIPADDR=$(grep ${HOSTNAME}-eth0 /etc/hosts | awk '{print $1}')
    if [ -n "$MICN" -a -n "$MICIPADDR" -a -n "$HOSTIPADDR" ]
    then
# generate /etc/mpss/default.conf files
      cat /etc/default.conf.template | sed "s/HOSTIPADDR/$HOSTIPADDR/g" > /tmp/default.conf
      delta=$(diff -u /tmp/default.conf /etc/mpss/default.conf | wc -l)
      if [ -n "$delta" -a $delta -gt 0 ]
      then
        cp /etc/mpss/default.conf /tmp/default.conf.bak
        cp /tmp/default.conf /etc/mpss/default.conf
      fi
# generate /etc/mpss/micN.conf files
      cat /etc/micN.conf.template | sed "s/NODE/$HOSTNAME/g" | sed "s/HOSTIPADDR/$HOSTIPADDR/g" | sed "s/MICN/$MICN/g" | sed "s/MICIPADDR/$MICIPADDR/g" > /tmp/${MICN}.conf
      delta=$(diff -u /tmp/${MICN}.conf /etc/mpss/${MICN}.conf | wc -l)
      if [ -n "$delta" -a $delta -gt 0 ]
      then
        cp /etc/mpss/${MICN}.conf /tmp/${MICN}.conf.bak
        cp /tmp/${MICN}.conf /etc/mpss/${MICN}.conf
      fi
    fi
  done
fi
