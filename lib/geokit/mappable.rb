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
  
      # Returns the centroid of a list of points, as a LatLng
      # Implementation of the algorithm described here:
      # http://www.jennessent.com/arcgis/shapes_poster.htm
      def centroid(locations, options = {})
        case locations.length
          when 1 then return locations.first
          when 2 then return midpoint_between(locations.first, locations.last, options)
          when 3 then return triangle_centroid(locations[0], locations[1], locations[2], options)
        end

        areas, total_area, centroids, clockwise = centroid_area_helper(locations, :all, options)

        # Take the area weighted sum of the centroids, with a negative weight
        # being applied to counter-clockwise trianges (which means its area is
        # entirely outside the polyogon)
        latitude = 0
        longitude = 0
        centroids.each_with_index do |centroid, i|
          weight = areas[i] / total_area
          weight *= -1 unless clockwise[i]
          latitude += centroid.lat * weight
          longitude += centroid.lng * weight
        end

        LatLng.new(latitude, longitude)
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
      def polygon(locations)
        locations = locations.dup
        start = locations.delete(northernmost(locations))

        headings = headings_between(start, locations)
        headings = headings.to_a
        headings.sort! { |a, b| a[1] <=> b[1] }

        [start] + headings.map { |v| v[0] }
      end
    
      def area(locations, options = {})
        if locations.length < 3
          # Not strictly true (see lunes), but tough luck
          0

        elsif locations.length == 3
          a, b, c = *locations

          # Use headings so that the sides determining angles are great circles
          ab, ac = heading_between(a, b), heading_between(a, c)
          ba, bc = heading_between(b, a), heading_between(b, c)
          ca, cb = heading_between(c, a), heading_between(c, b)
          angles = [(ab - ac).abs, (ba - bc).abs, (ca - cb).abs]

          total = angles.inject(0) do |sum, angle|
            angle = 360 - angle if angle > 180
            sum + angle / 180 * Math::PI
          end

          spherical_excess = total - Math::PI
          radius = units_sphere_multiplier(options[:units])

          spherical_excess * radius * radius

        else
          centroid_area_helper(locations, :area, options)
        end
      end

      def to_polar(x, y, z, options = {})
        radius = units_sphere_multiplier(options[:units])

        t = Math.acos(z / radius)
        lat = -1 * ((t * 180 / Math::PI) - 90)

        yx = y / x.to_f
        lng = if yx.nan?
          0
        else
          p = Math.atan(yx)
          l = p * 180 / Math::PI
          l += 180 if x < 0
          l
        end

        [lat, lng]
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

      def clockwise?(a, b, c)
        ab_heading = heading_between(a, b)
        ac_heading = heading_between(a, c)
        reverse_ab_heading = ab_heading == 180 ? 360 : ((ab_heading + 180) % 360)
        if ab_heading < reverse_ab_heading
          ac_heading > ab_heading && ac_heading < reverse_ab_heading
        else
          !(ac_heading < ab_heading && ac_heading > reverse_ab_heading)
        end
      end

      def triangle_centroid(a, b, c, options = {})
        radius = units_sphere_multiplier(options[:units])

        a = a.to_cartesian
        b = b.to_cartesian
        c = c.to_cartesian

        # Get the non-spherical centroid
        x = (a[0] + b[0] + c[0]) / 3
        y = (a[1] + b[1] + c[1]) / 3
        z = (a[2] + b[2] + c[2]) / 3

        # Normalize then multiply by the radius to move the point to the surface
        len = ((x ** 2) + (y ** 2) + (z ** 2)) ** 0.5
        x = x / len * radius
        y = y / len * radius
        z = z / len * radius

        LatLng.new(*to_polar(x, y, z))
      end

      # The calculations for the area of and the centroid of a set of points
      # larger than three are identical up to a point. That logic resides here,
      # with a parameter for which set of data is returned. Sort of messy, but
      # it keeps other things cleaner.
      def centroid_area_helper(locations, ret_value, options = {})
        # Can be :convex_hull or :polygon (be wary of the polygon method - see fixme below)
        polygon_method = options[:method] || :convex_hull
        polygon = send(polygon_method, locations)

        # This anchor point can actually be anywhere on the sphere but somewhere
        # in this implementation (or the implementation of the helpers), some
        # bug or presicion error is getting introduced over long distances that
        # results in the centroid being up to ~2 degrees off the mark. So, pick
        # an anchor inside or close to inside the polygon to minimize the error.
        # TODO: Find the bug described above
        # FIXME: This tries to be smart and works for convex hulls, but can pick
        #        a point that causes large inaccuracies for for non-convex
        #        polygons quite easily
        n = northernmost(polygon)
        i = polygon.index(n)
        mid = polygon[i - 1].midpoint_to(polygon[(i + 1) % polygon.length])
        anchor = n.endpoint(n.heading_to(mid), 1)

        # For each concecutive pair of verticies in the polygon, find the area
        # and centroid of each triangle formed by the pair and the anchor, and
        # the counter-clockwise or clockwise direction of the verticies.
        total_area, total_weighted_area = 0, 0
        areas, centroids, clockwise = [], [], []
        polygon << polygon.first
        polygon.each_cons(2) do |b, c|
          area = area([anchor, b, c])
          areas << area
          total_area += area

          cw = clockwise?(anchor, b, c)
          area *= -1 unless cw
          total_weighted_area += area
          clockwise << cw

          centroids << triangle_centroid(anchor, b, c, options)
        end

        if ret_value == :area
          total_weighted_area
        else
          [areas, total_area, centroids, clockwise]
        end
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

    # +X = 'towards', +Y = 'east', +Z = 'north', with latlng(0,0) being all 'towards'
    def to_cartesian(options = {})
      radius = self.class.send(:units_sphere_multiplier, options[:units])

      t = (90 - lat) * Math::PI / 180
      p = lng * Math::PI / 180

      sin_t = Math.sin(t)
      x = radius * sin_t * Math.cos(p)
      y = radius * sin_t * Math.sin(p)
      z = radius * Math.cos(t)

      [x, y, z]
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
