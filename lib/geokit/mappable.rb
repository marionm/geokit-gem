#require 'forwardable'

module Geokit     
  # Contains class and instance methods providing distance calcuation services.  This
  # module is meant to be mixed into classes containing lat and lng attributes where
  # distance calculation is desired.  
  # 
  # At present, two forms of distance calculations are provided:
  # 
  # * Pythagorean Theory (flat Earth) - which assumes the world is flat and loses accuracy over long distances.
  # * Haversine (sphere) - which is fairly accurate, but at a performance cost.
  # 
  # Distance units supported are :miles, :kms, and :nms.
  module Mappable
    PI_DIV_RAD = 0.0174
    KMS_PER_MILE = 1.609
    NMS_PER_MILE = 0.868976242
    EARTH_RADIUS_IN_MILES = 3963.19
    EARTH_RADIUS_IN_KMS = EARTH_RADIUS_IN_MILES * KMS_PER_MILE
    EARTH_RADIUS_IN_NMS = EARTH_RADIUS_IN_MILES * NMS_PER_MILE
    MILES_PER_LATITUDE_DEGREE = 69.1
    KMS_PER_LATITUDE_DEGREE = MILES_PER_LATITUDE_DEGREE * KMS_PER_MILE
    NMS_PER_LATITUDE_DEGREE = MILES_PER_LATITUDE_DEGREE * NMS_PER_MILE
    LATITUDE_DEGREES = EARTH_RADIUS_IN_MILES / MILES_PER_LATITUDE_DEGREE  
    
    # Mix below class methods into the includer.
    def self.included(receiver) # :nodoc:
      receiver.extend ClassMethods
    end   
    
    module ClassMethods #:nodoc:
      # Returns the distance between two points.  The from and to parameters are
      # required to have lat and lng attributes.  Valid options are:
      # :units - valid values are :miles, :kms, :nms (Geokit::default_units is the default)
      # :formula - valid values are :flat or :sphere (Geokit::default_formula is the default)
      def distance_between(from, to, options={})
        from=Geokit::LatLng.normalize(from)
        to=Geokit::LatLng.normalize(to)
        return 0.0 if from == to # fixes a "zero-distance" bug
        units = options[:units] || Geokit::default_units
        formula = options[:formula] || Geokit::default_formula
        case formula
        when :sphere
          begin
            units_sphere_multiplier(units) * 
                Math.acos( Math.sin(deg2rad(from.lat)) * Math.sin(deg2rad(to.lat)) + 
                Math.cos(deg2rad(from.lat)) * Math.cos(deg2rad(to.lat)) * 
                Math.cos(deg2rad(to.lng) - deg2rad(from.lng)))
          rescue Errno::EDOM
            0.0
          end
        when :flat
          Math.sqrt((units_per_latitude_degree(units)*(from.lat-to.lat))**2 + 
              (units_per_longitude_degree(from.lat, units)*(from.lng-to.lng))**2)
        end
      end

      # Returns heading in degrees (0 is north, 90 is east, 180 is south, etc)
      # from the first point to the second point. Typicaly, the instance methods will be used 
      # instead of this method.
      def heading_between(from,to)
        from=Geokit::LatLng.normalize(from)
        to=Geokit::LatLng.normalize(to)

        d_lng=deg2rad(to.lng-from.lng)
        from_lat=deg2rad(from.lat)
        to_lat=deg2rad(to.lat) 
        y=Math.sin(d_lng) * Math.cos(to_lat)
        x=Math.cos(from_lat)*Math.sin(to_lat)-Math.sin(from_lat)*Math.cos(to_lat)*Math.cos(d_lng)
        heading=to_heading(Math.atan2(y,x))
      end

      # Returns headings between a point and multiple others
      def headings_between(from, to)
        headings = {}
        to.each do |location|
          headings[location] = heading_between(from, location)
        end
        headings
      end
  
      # Given a start point, distance, and heading (in degrees), provides
      # an endpoint. Returns a LatLng instance. Typically, the instance method
      # will be used instead of this method.
      def endpoint(start,heading, distance, options={})
        units = options[:units] || Geokit::default_units
        radius = case units
          when :kms; EARTH_RADIUS_IN_KMS
          when :nms; EARTH_RADIUS_IN_NMS
          else EARTH_RADIUS_IN_MILES
        end
        start=Geokit::LatLng.normalize(start)        
        lat=deg2rad(start.lat)
        lng=deg2rad(start.lng)
        heading=deg2rad(heading)
        distance=distance.to_f
        
        end_lat=Math.asin(Math.sin(lat)*Math.cos(distance/radius) +
                          Math.cos(lat)*Math.sin(distance/radius)*Math.cos(heading))

        end_lng=lng+Math.atan2(Math.sin(heading)*Math.sin(distance/radius)*Math.cos(lat),
                               Math.cos(distance/radius)-Math.sin(lat)*Math.sin(end_lat))

        LatLng.new(rad2deg(end_lat),rad2deg(end_lng))
      end

      # Returns the midpoint, given two points. Returns a LatLng. 
      # Typically, the instance method will be used instead of this method.
      # Valid option:
      #   :units - valid values are :miles, :kms, or :nms (:miles is the default)
      def midpoint_between(from,to,options={})
        from=Geokit::LatLng.normalize(from)

        units = options[:units] || Geokit::default_units
        
        heading=from.heading_to(to)
        distance=from.distance_to(to,options)
        midpoint=from.endpoint(heading,distance/2,options)
      end
  
      # Geocodes a location using the multi geocoder.
      def geocode(location, options = {})
        res = Geocoders::MultiGeocoder.geocode(location, options)
        return res if res.success?
        raise Geokit::Geocoders::GeocodeError      
      end

      # Returns a clockwise-ordered subset of the locations that form a convex hull.
      # The algorithm is similar in concept to the polygon method's, but the heading
      # sorting logic is done for each added node rather than just the start point.
      # Won't really deal with map edge cases, so be careful.
      def convex_hull(locations)
        locations = locations.dup
        start = locations.delete(northernmost(locations))

        vertices = [start]
        current = nil
        last_heading = 0
        until locations.empty? do
          headings = headings_between(current || start, locations)
          locations << start unless current

          # Force headings to be sorted in 'clockwise-closeness' order to the
          # heading between the previous point and the current. This prevents
          # any inner point from getting picked inappropriately.
          headings.each do |location, heading|
            headings[location] += 360 if heading < last_heading
          end

          headings = headings.to_a
          headings.sort! { |a, b| a[1] <=> b[1] }

          current = headings.first[0]
          last_heading = headings.first[1]

          break if current == start

          vertices << current
          locations.delete(current)
        end

        vertices
      end

      # Returns a clockwise ordering of the locations that will form a polygon.
      # Won't really deal with map edge cases, so be careful.
      def polygon(locations)
        locations = locations.dup
        start = locations.delete(northernmost(locations))

        headings = headings_between(start, locations)
        headings = headings.to_a
        headings.sort! { |a, b| a[1] <=> b[1] }

        [start] + headings.map { |v| v[0] }
      end
    
      protected
    
      def deg2rad(degrees)
        degrees.to_f / 180.0 * Math::PI
      end
      
      def rad2deg(rad)
        rad.to_f * 180.0 / Math::PI 
      end
      
      def to_heading(rad)
        (rad2deg(rad)+360)%360
      end

      # Returns the multiplier used to obtain the correct distance units.
      def units_sphere_multiplier(units)
        case units
          when :kms; EARTH_RADIUS_IN_KMS
          when :nms; EARTH_RADIUS_IN_NMS
          else EARTH_RADIUS_IN_MILES
        end
      end

      # Returns the number of units per latitude degree.
      def units_per_latitude_degree(units)
        case units
          when :kms; KMS_PER_LATITUDE_DEGREE
          when :nms; NMS_PER_LATITUDE_DEGREE
          else MILES_PER_LATITUDE_DEGREE
        end
      end
    
      # Returns the number units per longitude degree.
      def units_per_longitude_degree(lat, units)
        miles_per_longitude_degree = (LATITUDE_DEGREES * Math.cos(lat * PI_DIV_RAD)).abs
        case units
          when :kms; miles_per_longitude_degree * KMS_PER_MILE
          when :nms; miles_per_longitude_degree * NMS_PER_MILE
          else miles_per_longitude_degree
        end
      end  

      # Return the northern-most point from the list. Ties are won by a lower array index.
      def northernmost(locations)
        max = locations.first
        locations[1..-1].each do |location|
          max = location if location.lat > max.lat
        end
        max
      end
    end
  
    # -----------------------------------------------------------------------------------------------
    # Instance methods below here
    # -----------------------------------------------------------------------------------------------
  
    # Extracts a LatLng instance. Use with models that are acts_as_mappable
    def to_lat_lng
      return self if instance_of?(Geokit::LatLng) || instance_of?(Geokit::GeoLoc)
      return LatLng.new(send(self.class.lat_column_name),send(self.class.lng_column_name)) if self.class.respond_to?(:acts_as_mappable)
      nil
    end

    # Returns the distance from another point.  The other point parameter is
    # required to have lat and lng attributes.  Valid options are:
    # :units - valid values are :miles, :kms, :or :nms (:miles is the default)
    # :formula - valid values are :flat or :sphere (:sphere is the default)
    def distance_to(other, options={})
      self.class.distance_between(self, other, options)
    end  
    alias distance_from distance_to

    # Returns heading in degrees (0 is north, 90 is east, 180 is south, etc)
    # to the given point. The given point can be a LatLng or a string to be Geocoded 
    def heading_to(other)
      self.class.heading_between(self,other)
    end

    # Returns heading in degrees (0 is north, 90 is east, 180 is south, etc)
    # FROM the given point. The given point can be a LatLng or a string to be Geocoded 
    def heading_from(other)
      self.class.heading_between(other,self)
    end
 
    # Returns the endpoint, given a heading (in degrees) and distance.  
    # Valid option:
    # :units - valid values are :miles, :kms, or :nms (:miles is the default)
    def endpoint(heading,distance,options={})
      self.class.endpoint(self,heading,distance,options)  
    end

    # Returns the midpoint, given another point on the map.  
    # Valid option:
    # :units - valid values are :miles, :kms, or :nms (:miles is the default)    
    def midpoint_to(other, options={})
      self.class.midpoint_between(self,other,options)
    end
    
  end
end
