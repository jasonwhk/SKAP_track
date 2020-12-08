import threading
import time
import random
import tornado
import signal
import json
import time
from astropy.time import Time
from katcp import AsyncDeviceServer, Message, DeviceServer, Sensor, ProtocolFlags, AsyncReply
from katcp.kattypes import (Str, Float, Timestamp, Discrete, request, return_reply)
from katcp.resource_client import KATCPClientResource

import sys

server_host = "0.0.0.0"

print(sys.argv[1])
server_port = int(sys.argv[1])


class MyServer(AsyncDeviceServer):

    VERSION_INFO = ("example-api", 1, 0)
    BUILD_INFO = ("example-implementation", 0, 1, "")

    # Optionally set the KATCP protocol version and features. Defaults to
    # the latest implemented version of KATCP, with all supported optional
    # featuresthat's all of the receivers
    PROTOCOL_INFO = ProtocolFlags(5, 0, set([
        ProtocolFlags.MULTI_CLIENT,
        ProtocolFlags.MESSAGE_IDS,
    ]))

    FRUIT = [
        "apple", "banana", "pear", "kiwi",
    ]

    def setup_sensors(self):
        """Setup some server sensors."""
	
        self._add_result = Sensor.float("add.result",
                                        "Last ?add result.", "", [-10000, 10000])

        self._time_result = Sensor.timestamp("time.result",
                                             "Last ?time result.", "")

        self._eval_result = Sensor.string("eval.result",
                                          "Last ?eval result.", "")

        self._fruit_result = Sensor.discrete("fruit.result",
                                             "Last ?pick-fruit result.", "", self.FRUIT)
        self._device_armed = Sensor.boolean(
            "device-armed",
            description="Is the CAM server armed?",
            initial_status=Sensor.NOMINAL,
            default=True)
        self._bandwidth = Sensor.float("bandwidth", default=300)
        self._sourcename = Sensor.string("sourcename", default="none")
        self._source_ra = Sensor.string("source_RA", default=0)
        self._source_dec = Sensor.string("source_DEC", default=0)
        self._exposure_time = Sensor.float("EXP_time", default=0)

        self.add_sensor(self._sourcename)
        self.add_sensor(self._source_ra)
        self.add_sensor(self._source_dec)
        self.add_sensor(self._exposure_time)

        self.add_sensor(self._bandwidth)
        self.add_sensor(self._device_armed)
        self.add_sensor(self._add_result)
        self.add_sensor(self._time_result)
        self.add_sensor(self._eval_result)
        self.add_sensor(self._fruit_result)

        self._systemp_result = Sensor.float("add.result",
                                            "Last ?add result.", "", [-10000, 10000])
        self.add_sensor(self._systemp_result)
	
	##self._bandwidth = Sensor.float("bandwidth", default=300)
	#self.add_sensor(self._bandwidth)

    @request()
    @return_reply(Str())
    def request_bandwidth(self, req, bw):
        """Return the Bandwidth"""
        #req.inform("checking armed status", self._device_armed.value())
        req.reply("ok", bw)
        raise AsyncReply
    @request()
    @return_reply(Str())
    def request_status_armed(self, req):
        """Return the state of the Armed/Disarmed"""
        req.inform("checking armed status", self._device_armed.value())
        req.reply("ok", self._device_armed.value())
        raise AsyncReply

    @request(Float())
    @return_reply()
    def request_long_action(self, req, t):
        """submit a long action command for testing using coroutine"""
        @tornado.gen.coroutine
        def wait():
            yield tornado.gen.sleep(t)
            req.reply("slept for", t, "second")
        self.ioloop.add_callback(wait)
        raise AsyncReply

    @request(Float(), Float())
    @return_reply(Str())
    def request_radec(self, req, ra, dec):
        """testing to read in the RA DEC fomr a client"""
        # test=ra+dec
        self.ra = ra
        self.dec = dec
        return ("ok", "%f %f" % (self.ra, self.dec))

    @request(Float(), Float())
    @return_reply(Float())
    def request_add(self, req, x, y):
        """Add two numbers"""
        r = x + y
        self._add_result.set_value(r)
        return ("ok", r)

    @request()
    @return_reply(Str())
    def request_arm(self, req):
        """Arm the controller"""
        @tornado.gen.coroutine
        def start_controller():
            req.inform("processing", "command processing")
            try:
                yield tornado.gen.sleep(10)
            except Exception as error:
                req.reply("fail", "Unknown error: {0}".format(str(error)))
            else:
                req.reply("ok", "effcam armed")
                self._device_armed.set_value(True)
        if self._device_armed.value():
            return ("fail", "Effcam is already armed")
        self.ioloop.add_callback(start_controller)
        raise AsyncReply

    @request()
    @return_reply(Str())
    def request_disarm(self, req):
        """disarm the controller"""
        @tornado.gen.coroutine
        # @coroutine
        def stop_controller():
            req.inform("processing", "processing command")
            try:
                yield tornado.gen.sleep(10)
                # yield self._controller.stop()
            except Exception as error:
                req.reply("fail", "Unknown error: {0}".format(str(error)))
            else:
                req.reply("ok", "effcam disarmed")
                self._device_armed.set_value(False)
        if self._device_armed.value() == False:
            return ("fail", "Effcam is already disarmed")
        self.ioloop.add_callback(stop_controller)
        raise AsyncReply

    @request()
    @return_reply(Str())
    def request_status_temp(self, req):
        """Return the current temp"""
        #r = time.time()
        t = "36"
        # self._time_result.set_value(r)
        return ("ok", t)

    @request()
    @return_reply(Timestamp())
    def request_status_time(self, req):
        """Return the current time in seconds since the Unix Epoch."""
        req.inform("processing", "processing command")
        r = time.time()
        # self._time_result.set_value(r)
        req.reply("ok", r)
        raise AsyncReply
        # return ("ok", r)

    @request()
    @return_reply(Timestamp(), Str())
    def request_status_time_and_temp(self, req):
        """Return the current time in seconds since the Unix Epoch."""
        req.inform("processing", "processing command")
        r = time.time()
        # self._time_result.set_value(r)
        t = "36"
        req.reply("ok", r, t)
        raise AsyncReply

    @request(Str())
    @return_reply()
    def request_configure(self, req, config):
        """Return ok."""
	print "{} received configuration {}".format(Time.now(),config)
        self.config = config
	time.sleep(1)
	req.reply("ok",)
        raise AsyncReply

    @request(Str())
    @return_reply()
    def request_provision(self, req, config):
        """Return ok."""
        print "{} received provision {}".format(Time.now(),config)
        self.config = config
        time.sleep(1)
        req.reply("ok",)
        raise AsyncReply

    @request(Str())
    @return_reply()
    def request_measurement_prepare(self, req, config):
        """Return ok."""
        print "{} received measurement prepare {}".format(Time.now(),config)
        self.config = config
        time.sleep(1)
        req.reply("ok",)
        raise AsyncReply

    @request(Str())
    @return_reply()
    def request_configure(self, req, config):
        """Return ok."""
        print "{} received configuration {}".format(Time.now(),config)
        self.config = config
        time.sleep(1)
        req.reply("ok",)
        raise AsyncReply

    @request()
    @return_reply(Str())
    def request_status_config(self, req):
        """Return ok."""
        req.reply("ok", "{}".format(self.config))
        raise AsyncReply

    @request()
    @return_reply()
    def request_capture_start(self, req):
        """Return ok."""
        print "{} received capture start request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply

    @request()
    @return_reply()
    def request_capture_stop(self, req):
        """Return ok."""
        print "{} received capture stop request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply

    @request()
    @return_reply()
    def request_measurement_start(self, req):
        """Return ok."""
        print "{} received measurement start request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply

    @request()
    @return_reply()
    def request_measurement_stop(self, req):
        """Return ok."""
        print "{} received measurement stop request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply

    @request()
    @return_reply()
    def request_deconfigure(self, req):
        """Return ok."""
        print "{} received deconfigure request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply

    @request()
    @return_reply()
    def request_deprovision(self, req):
        """Return ok."""
        print "{} received deprovision request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply()

    @return_reply()
    def request_start(self, req):
        """Return ok."""
	print "{} received start request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply

    @request()
    @return_reply()
    def request_stop(self, req):
        """Return ok."""
	print "{} received stop request on port :{}".format(Time.now(), server_port)
        req.reply("ok")
        raise AsyncReply


@tornado.gen.coroutine
def on_shutdown(ioloop, server):
    print('Shutting down')
    yield server.stop()
    ioloop.stop()


if __name__ == "__main__":
    ioloop = tornado.ioloop.IOLoop.current()
    #logging.getLogger().addHandler(logging.NullHandler())
    #logger = logging.getLogger('dummy_backend')
    #coloredlogs.install(
    #    fmt="[ %(levelname)s - %(asctime)s - %(name)s - %(filename)s:%(lineno)s] %(message)s",
    #    level="INFO",
    #    logger=logger)
    #logging.getLogger('dummy_backend').setLevel('INFO')
    server = MyServer(server_host, server_port)
    signal.signal(signal.SIGINT, lambda sig, frame: ioloop.add_callback_from_signal(
        on_shutdown, ioloop, server))
    ioloop.add_callback(server.start)
    print "{} Starting dummy backend instance".format(Time.now())
    ioloop.start()
